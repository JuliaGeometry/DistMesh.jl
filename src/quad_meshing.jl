################################################################################
# quad_meshing.jl  –  Triangle-to-quad conversion for DistMesh.jl
#
# Pipeline:
#   Stage 1: Topology-driven tri-to-quad matching (Blossom perfect matching)
#   Stage 2: Boundary triangle collapse (union-find node merging)
#   Stage 3: Greedy quad topology improvement (diagonal collapse + edge flip)
#   Stage 4: Optional Catmull-Clark refinement (if triangles survive stage 2)
################################################################################

################################################################################
# Internal topology cache
################################################################################

struct _QuadTopology
    is_corner      :: Vector{Bool}
    is_boundary    :: Vector{Bool}
    nb             :: Matrix{Tuple{Int32,Int32}}
    element_depths :: Vector{Int}
    node_depths    :: Vector{Int}
    adj            :: Vector{Vector{Int}}
    degree         :: Vector{Int}
    ideal_degree   :: Vector{Int}
    min_degree     :: Vector{Int}
end

# Match pfix corner coordinates to nearest mesh nodes
function _corner_node_ids(p::Vector{<:SVector}, pfix)
    isempty(pfix) && return Int[]
    corners = Int[]
    for pc in pfix
        pcsv = SVector{2,Float64}(Float64(pc[1]), Float64(pc[2]))
        idx = argmin(norm(pi - pcsv) for pi in p)
        push!(corners, idx)
    end
    return unique!(sort!(corners))
end

# BFS element depths from boundary elements (depth 0 = on boundary)
function _element_depths(nb::Matrix{Tuple{Int32,Int32}})
    nt = size(nb, 2)
    nf = size(nb, 1)
    D = fill(-1, nt)
    queue = Int[]
    for it in 1:nt
        for jf in 1:nf
            if nb[jf, it][1] == 0
                D[it] = 0
                push!(queue, it)
                break  # break inner loop only; continue outer loop
            end
        end
    end
    for current in queue
        for jf in 1:nf
            jel = Int(nb[jf, current][1])
            if jel > 0 && D[jel] == -1
                D[jel] = D[current] + 1
                push!(queue, jel)
            end
        end
    end
    return D
end

# BFS node depths from boundary nodes (depth 0 = on boundary)
function _node_depths(is_boundary::Vector{Bool}, adj::Vector{Vector{Int}})
    D = fill(-1, length(is_boundary))
    queue = findall(is_boundary)
    D[queue] .= 0
    for current in queue
        for nb in adj[current]
            if D[nb] == -1
                D[nb] = D[current] + 1
                push!(queue, nb)
            end
        end
    end
    return D
end

# Compute ideal and minimum degrees for each node.
# Interior: ideal=4, min=3
# Smooth boundary (not in pfix): ideal=3, min=2
# Corner (in pfix): angle-based ideal and min
function _target_degrees(msh::DMesh, is_boundary::Vector{Bool}, is_corner::Vector{Bool};
                         angle_tol=10.0)
    np = length(msh.p)
    ideal_d = fill(4, np)
    min_d   = fill(3, np)

    # Smooth boundary nodes
    for i in 1:np
        if is_boundary[i] && !is_corner[i]
            ideal_d[i] = 3
            min_d[i]   = 3
        end
    end

    # Corner nodes: compute interior angle from ordered boundary segments
    if any(is_corner)
        segs = segcollect(boundary_faces(msh))
        for seg in segs
            N = length(seg)
            for i in 1:N
                curr = seg[i]
                is_corner[curr] || continue
                prev = seg[mod1(i-1, N)]
                nxt  = seg[mod1(i+1, N)]
                v_in  = msh.p[curr] - msh.p[prev]
                v_out = msh.p[nxt]  - msh.p[curr]
                theta = 180.0 - rad2deg(atan(v_in[1]*v_out[2] - v_in[2]*v_out[1],
                                             v_in[1]*v_out[1] + v_in[2]*v_out[2]))
                ideal_d[curr] = max(round(Int, theta / 90.0) + 1, 2)
                max_angle = 180.0 - angle_tol
                if theta < max_angle
                    min_d[curr] = 2
                elseif theta < 2.0 * max_angle
                    min_d[curr] = 3
                else
                    min_d[curr] = 4
                end
            end
        end
    end

    return ideal_d, min_d
end

function _QuadTopology(msh::DMesh, pfix)
    bnodes = boundary_nodes(msh)
    np = length(msh.p)
    is_boundary = falses(np)
    is_boundary[bnodes] .= true

    corner_ids = _corner_node_ids(msh.p, pfix)
    is_corner = falses(np)
    is_corner[corner_ids] .= true

    nb = element_face_neighbors(msh)

    adj = [Int.(a) for a in node_adjacency(msh)]
    deg = node_degrees(msh)
    ideal_d, min_d = _target_degrees(msh, Vector{Bool}(is_boundary), Vector{Bool}(is_corner))

    return _QuadTopology(
        Vector{Bool}(is_corner),
        Vector{Bool}(is_boundary),
        nb,
        _element_depths(nb),
        _node_depths(Vector{Bool}(is_boundary), adj),
        adj,
        deg,
        ideal_d,
        min_d,
    )
end

################################################################################
# BlossomV perfect matching wrapper
################################################################################

# g: Vector of (u, v, weight) 0-indexed element pairs
# N: number of elements
# Returns 1-indexed matching; out[i]=j means i matched with j; 0 = unmatched
function _blossomv_match(g, N)
    matching = Matching(Int, 2N, 2 * length(g) + N)
    for (u, v, w) in g
        add_edge(matching, u, v, -w)           # real edge (negated for min-cost)
        add_edge(matching, u + N, v + N, 0)    # dummy mirror
    end
    for i in 0:(N-1)
        add_edge(matching, i, i + N, 0)        # escape: allow unmatched
    end
    solve(matching)
    out = get_match.([matching], 0:N-1) .+ 1
    out[out .> N] .= 0
    return out
end

################################################################################
# Stage 1: Triangle-to-quad matching
################################################################################

function _dual_mesh_edges_weighted(msh::DMesh, topo::_QuadTopology, nb::Matrix)
    emap = edgemap(Simplex{2}())
    nt   = length(msh.t)
    nf   = 3

    d_curr      = topo.degree
    d_ideal     = topo.ideal_degree
    d_min       = topo.min_degree
    D           = topo.element_depths
    is_boundary = topo.is_boundary
    is_corner   = topo.is_corner

    w_base              = 0.0
    gamma               = 1000.0
    beta_scale          = 1.0
    boundary_strictness = 100.0

    g = Tuple{Int,Int,Int}[]

    for it in 1:nt, ie in 1:nf
        jt = Int(nb[ie, it][1])
        jt > it || continue  # skip boundary (jt=0) and already-processed pairs

        weight = w_base + gamma * (D[it] + D[jt])

        u = Int(msh.t[it][emap[ie][1]])
        v = Int(msh.t[it][emap[ie][2]])

        gain_u = 2 * (d_curr[u] - d_ideal[u])
        gain_v = 2 * (d_curr[v] - d_ideal[v])
        mult_u = is_boundary[u] ? boundary_strictness : 1.0
        mult_v = is_boundary[v] ? boundary_strictness : 1.0

        # Skip if merging would drop a node below its minimum degree
        (d_curr[u] - 1 < d_min[u] || d_curr[v] - 1 < d_min[v]) && continue

        weight += beta_scale * (mult_u * gain_u + mult_v * gain_v)
        push!(g, (it - 1, jt - 1, round(Int, weight)))
    end

    return g
end

"""
    match_tri2quad(msh::DMesh, pfix) -> (q, t, t_global_idx)

Stage 1: Match pairs of adjacent triangles into quads via Blossom perfect matching
on the weighted dual graph. Returns quad connectivity, unmatched triangle
connectivity, and the global indices of unmatched triangles.
"""
function match_tri2quad(msh::DMesh, topo)
    nb   = topo.nb
    nt   = length(msh.t)

    g   = _dual_mesh_edges_weighted(msh, topo, nb)
    out = _blossomv_match(g, nt)

    emap = edgemap(Simplex{2}())
    q = Index4[]

    for it in 1:nt
        jt = out[it]
        jt > it || continue  # process each matched pair once

        # Which local face of it borders jt, and vice versa
        j = findfirst(ie -> Int(nb[ie, it][1]) == jt, 1:3)
        k = Int(nb[j, it][2])

        # Quad nodes: "next edge" of it, then "next edge" of jt
        next_j = j % 3 + 1
        next_k = k % 3 + 1
        n1 = Int(msh.t[it][emap[next_j][1]])
        n2 = Int(msh.t[it][emap[next_j][2]])
        n3 = Int(msh.t[jt][emap[next_k][1]])
        n4 = Int(msh.t[jt][emap[next_k][2]])
        push!(q, Index4(n1, n2, n3, n4))
    end

    t_global_idx = findall(out .== 0)
    t = [msh.t[i] for i in t_global_idx]
    return q, t, t_global_idx
end

################################################################################
# Stage 2: Boundary triangle collapse
################################################################################

# Node degrees for a mixed quad+triangle mesh, counting each shared edge once.
function _mixed_node_degrees(np::Int, q::Vector{Index4}, t::Vector{<:SVector{3}})
    edges = Set{Tuple{Int,Int}}()
    for quad in q
        a, b, c, d = Int(quad[1]), Int(quad[2]), Int(quad[3]), Int(quad[4])
        push!(edges, minmax(a,b)); push!(edges, minmax(b,c))
        push!(edges, minmax(c,d)); push!(edges, minmax(d,a))
    end
    for tri in t
        a, b, c = Int(tri[1]), Int(tri[2]), Int(tri[3])
        push!(edges, minmax(a,b)); push!(edges, minmax(b,c)); push!(edges, minmax(c,a))
    end
    deg = zeros(Int, np)
    for (a, b) in edges; deg[a] += 1; deg[b] += 1; end
    return deg
end

function _collapse_boundary_triangles(msh::DMesh, topo,
                                      q::Vector{Index4},
                                      t_unmatched::Vector{<:SVector{3}},
                                      t_global_idx::Vector{Int})
    np         = length(msh.p)
    emap       = edgemap(Simplex{2}())
    nb         = topo.nb
    is_corner  = topo.is_corner
    min_degree = topo.min_degree
    degree     = _mixed_node_degrees(np, q, t_unmatched)

    # Union-find with corner pinning:
    # corner_root[i] != 0 means root i represents the pfix corner at that index.
    parent      = collect(1:np)
    corner_root = [is_corner[i] ? i : 0 for i in 1:np]

    find_root(i) = begin
        while parent[i] != i; i = parent[i]; end
        i
    end

    function union!(i, j)
        ri, rj = find_root(i), find_root(j)
        ri == rj && return
        (corner_root[ri] != 0 && corner_root[rj] != 0) && return  # both corners: skip
        if corner_root[rj] != 0
            ri, rj = rj, ri  # ensure corner becomes root
        end
        parent[rj] = ri
        corner_root[ri] = max(corner_root[ri], corner_root[rj])
    end

    # For each unmatched triangle, collapse its first valid boundary edge:
    # valid = the opposite node retains degree >= its minimum after losing this triangle.
    # Degrees are updated incrementally so later decisions see accurate values.
    for gt in t_global_idx
        for ie in 1:3
            Int(nb[ie, gt][1]) == 0 || continue
            a        = Int(msh.t[gt][emap[ie][1]])
            b        = Int(msh.t[gt][emap[ie][2]])
            opp      = Int(msh.t[gt][ie])   # local node opposite to face ie in Simplex{2}
            opp_root = find_root(opp)
            degree[opp_root] - 1 >= min_degree[opp] || continue
            ra = find_root(a)
            rb = find_root(b)
            union!(a, b)
            r_merged = find_root(a)
            r_other  = r_merged == ra ? rb : ra
            degree[r_merged] = degree[ra] + degree[rb] - 2  # gt removed; b's elements absorbed
            degree[r_other]  = 0                            # no longer a root
            degree[opp_root] -= 1                           # gt removed from opp's neighbourhood
            break
        end
    end

    # Path-compress, then assign contiguous new IDs to the surviving roots.
    for i in 1:np; parent[i] = find_root(i); end
    unique_roots = unique(parent)
    root_to_new  = Dict(r => id for (id, r) in enumerate(unique_roots))
    old_to_new   = [root_to_new[parent[i]] for i in 1:np]

    # New node positions: corners snap to their pfix coordinate; others are
    # the centroid of merged members (single O(n) pass, no per-root findall).
    p_sum = fill(zero(eltype(msh.p)), length(unique_roots))
    p_cnt = zeros(Int, length(unique_roots))
    for i in 1:np
        nid = old_to_new[i]
        p_sum[nid] += msh.p[i]
        p_cnt[nid] += 1
    end
    new_p = [corner_root[r] != 0 ? msh.p[corner_root[r]] : p_sum[id] / p_cnt[id]
             for (id, r) in enumerate(unique_roots)]

    # Remap connectivity; drop degenerate (repeated-node) elements.
    new_q = filter(rq -> length(unique(rq)) == 4,
                   [Index4(old_to_new[quad[1]], old_to_new[quad[2]],
                           old_to_new[quad[3]], old_to_new[quad[4]]) for quad in q])
    new_t = filter(rt -> length(unique(rt)) == 3,
                   [Index3(old_to_new[tri[1]], old_to_new[tri[2]],
                           old_to_new[tri[3]]) for tri in t_unmatched])

    return new_p, new_q, new_t
end

################################################################################
# Stage 3: Catmull-Clark refinement (mixed tri/quad → all-quad)
################################################################################

function _catmull_clark_refine(p::Vector{Point2d}, q::Vector{Index4}, t::Vector{Index3})
    qemap = edgemap(Block{2}())
    temap = edgemap(Simplex{2}())
    nq, nt = length(q), length(t)

    # Deduplicated edge midpoints
    emid_dict = Dict{Tuple{Int,Int}, Int}()
    emids = Point2d[]

    function get_or_add_emid(a::Int, b::Int)
        key = a < b ? (a, b) : (b, a)
        get!(emid_dict, key) do
            push!(emids, (p[a] + p[b]) / 2)
            length(emids)
        end
    end

    for tri in t,  e_local in temap; get_or_add_emid(Int(tri[e_local[1]]), Int(tri[e_local[2]])); end
    for quad in q, e_local in qemap; get_or_add_emid(Int(quad[e_local[1]]), Int(quad[e_local[2]])); end

    np = length(p)
    ne = length(emids)

    # Element centroids (triangles first, then quads)
    t_cents = [sum(p[Int(i)] for i in tri)  / 3  for tri  in t]
    q_cents = [sum(p[Int(i)] for i in quad) / 4  for quad in q]

    new_quads = Index4[]

    # Each triangle → 3 quads
    # temap edges: [1]=(2,3), [2]=(3,1), [3]=(1,2)
    for (it, tri) in enumerate(t)
        em = [get_or_add_emid(Int(tri[e_local[1]]), Int(tri[e_local[2]])) for e_local in temap]
        mid = np + ne + it
        push!(new_quads, Index4(tri[1], np+em[3], mid, np+em[2]))
        push!(new_quads, Index4(tri[2], np+em[1], mid, np+em[3]))
        push!(new_quads, Index4(tri[3], np+em[2], mid, np+em[1]))
    end

    # Each quad → 4 quads
    # qemap edges: [1]=(1,2), [2]=(2,3), [3]=(3,4), [4]=(4,1)
    for (iq, quad) in enumerate(q)
        em = [get_or_add_emid(Int(quad[e_local[1]]), Int(quad[e_local[2]])) for e_local in qemap]
        mid = np + ne + nt + iq
        push!(new_quads, Index4(quad[1], np+em[1], mid, np+em[4]))
        push!(new_quads, Index4(quad[2], np+em[2], mid, np+em[1]))
        push!(new_quads, Index4(quad[3], np+em[3], mid, np+em[2]))
        push!(new_quads, Index4(quad[4], np+em[4], mid, np+em[3]))
    end

    new_p = vcat(p, emids, t_cents, q_cents)
    return DMesh(new_p, new_quads)
end

################################################################################
# Stage 4: Greedy quad topology operations (collapse + edge flip)
################################################################################

function _delta_E(a, b, c, d, deg, ideal_deg, D)
    W(n) = 100 - D[n]
    E(n, dnew) = W(n) * (dnew - ideal_deg[n])^2
    E_old = E(a, deg[a]) + E(b, deg[b]) + E(c, deg[c]) + E(d, deg[d])
    E_new = E(a, deg[a] + deg[c] - 2) + E(b, deg[b] - 1) + E(d, deg[d] - 1)
    return E_new - E_old
end

# Energy change for rotating the shared edge (n0,n3) → (n_lose1, n_lose2 lose 1 each;
# n_gain1, n_gain2 gain 1 each).
function _delta_E_flip(n_lose1, n_lose2, n_gain1, n_gain2, deg, ideal_deg, D)
    W(n) = 100 - D[n]
    E(n, dnew) = W(n) * (dnew - ideal_deg[n])^2
    nodes = (n_lose1, n_lose2, n_gain1, n_gain2)
    E_old = sum(E(n, deg[n]) for n in nodes)
    E_new = E(n_lose1, deg[n_lose1]-1) + E(n_lose2, deg[n_lose2]-1) +
            E(n_gain1, deg[n_gain1]+1) + E(n_gain2, deg[n_gain2]+1)
    return E_new - E_old
end

function _collapse_safe(a, c, is_boundary, is_corner, adj)
    is_boundary[a] && is_boundary[c]                       && return false
    (is_corner[a] || is_corner[c])                         && return false
    shared = intersect(adj[a], adj[c])
    length(shared) != 2                                    && return false
    (is_boundary[a] || is_boundary[c]) &&
        (is_boundary[shared[1]] || is_boundary[shared[2]]) && return false
    return true
end

# Safety check for a quad edge flip: the two shared-edge endpoints lose degree,
# so they must not go below min_degree. No node is merged, so corners are allowed
# to participate as long as they don't go below their minimum.
function _flip_safe(n0, n3, is_boundary, is_corner, min_degree, deg)
    deg[n0] - 1 < min_degree[n0] && return false
    deg[n3] - 1 < min_degree[n3] && return false
    return true
end

function _greedy_quad_improve(msh::DMesh, pfix; verbose=false)
    topo = _QuadTopology(msh, pfix)
    np = length(msh.p)

    p           = copy(msh.p)
    quads       = [[Int(n) for n in q] for q in msh.t]  # mutable Vector{Vector{Int}}
    is_dead     = falses(length(quads))
    deg         = copy(topo.degree)
    ideal_deg   = topo.ideal_degree
    min_degree  = topo.min_degree
    is_boundary = topo.is_boundary
    is_corner   = topo.is_corner
    live        = collect(1:length(quads))
    tmp         = msh
    adj         = [Int.(a) for a in node_adjacency(tmp)]
    D           = _node_depths(is_boundary, adj)
    nb_live     = element_face_neighbors(tmp)  # (4, nq) matrix, indexed by live position
    final_node  = collect(1:np)  # tracks which node each node was merged into

    for pass in 1:100
        # ── Candidate: best collapse ──────────────────────────────────────────
        best_dE   = 0.0
        best_op   = :none   # :collapse or :flip
        # collapse fields
        best_a = best_c = best_b = best_d = 0
        # flip fields: n0..n5 in canonical ring order; original quads indices i1, i2
        best_flip_n  = (0,0,0,0,0,0)
        best_flip_i1 = best_flip_i2 = 0
        best_flip_steps = 0

        for i in 1:length(quads)
            is_dead[i] && continue
            ns = quads[i]
            # Two diagonal pairs per quad
            for (ai, bi, ci, di) in ((1,2,3,4), (2,3,4,1))
                a, b, c, d = ns[ai], ns[bi], ns[ci], ns[di]
                _collapse_safe(a, c, is_boundary, is_corner, adj) || continue
                dE = _delta_E(a, b, c, d, deg, ideal_deg, D)
                if dE < best_dE
                    best_dE = dE
                    best_op = :collapse
                    best_a  = is_boundary[a] ? a : (is_boundary[c] ? c : a)
                    best_c  = (best_a == a) ? c : a
                    best_b, best_d = b, d
                end
            end
        end

        # ── Candidate: best edge flip ─────────────────────────────────────────
        # nb_live is indexed by live-quad position; translate to original quads[] index
        # via the `live` array: live[li] = original index i.
        for (li, i) in enumerate(live)
            ns = quads[i]
            for k in 1:4
                lj = Int(nb_live[k, li][1])
                lj == 0 && continue              # boundary
                j = live[lj]
                j <= i && continue               # process each pair once
                kj = Int(nb_live[k, li][2])

                ms = quads[j]

                # Build the canonical 6-node ring for this two-quad patch.
                # Quad i winding: [ea, eb, ec, ed], where ea-eb is the shared face edge.
                # Quad j shares the same edge (ea-eb or eb-ea).
                # Ring order around the combined patch: ea, ed, ec, eb, fc, fd
                # (n0=ea, n1=ed, n2=ec, n3=eb, n4=fc, n5=fd)
                # where fc, fd are the two private nodes of quad j, ordered so that
                # in quad j the winding goes ea → fd → fc → eb (i.e., fc is adjacent to eb).
                #
                # Flips (no node merging; two nodes gain/lose 1 degree each):
                #   1-step: new shared = n1-n4 = ed-fc  →  [ea,ed,fc,fd] + [ed,ec,eb,fc]
                #   2-step: new shared = n2-n5 = ec-fd  →  [ea,ed,ec,fd] + [ec,eb,fc,fd]
                ea = ns[k]
                eb = ns[mod1(k + 1, 4)]
                ec = ns[mod1(k + 2, 4)]
                ed = ns[mod1(k + 3, 4)]

                fa = ms[kj]  # = ea or eb
                if fa == ea
                    # quad j winding: [ea, eb, fc, fd] → fc adj to eb, fd adj to ea
                    fc_node = ms[mod1(kj + 2, 4)]
                    fd_node = ms[mod1(kj + 3, 4)]
                else
                    # quad j winding: [eb, ea, X, Y] → X adj to ea = fd, Y adj to eb = fc
                    fc_node = ms[mod1(kj + 3, 4)]
                    fd_node = ms[mod1(kj + 2, 4)]
                end

                n = (ea, ed, ec, eb, fc_node, fd_node)  # n[1..6] = n0..n5

                for (steps, lose1, lose2, gain1, gain2) in (
                        (1, n[1], n[4], n[2], n[5]),   # 1-step: n0,n3 lose; n1,n4 gain
                        (2, n[1], n[4], n[3], n[6]))   # 2-step: n0,n3 lose; n2,n5 gain
                    _flip_safe(lose1, lose2, is_boundary, is_corner, min_degree, deg) || continue
                    dE = _delta_E_flip(lose1, lose2, gain1, gain2, deg, ideal_deg, D)
                    if dE < best_dE
                        best_dE         = dE
                        best_op         = :flip
                        best_flip_n     = n
                        best_flip_i1    = i  # original quads[] index
                        best_flip_i2    = j  # original quads[] index
                        best_flip_steps = steps
                    end
                end
            end
        end

        best_dE >= -0.01 && (verbose && println("  Quad improve converged in $pass passes."); break)

        if best_op == :collapse
            a, c = best_a, best_c
            if !is_boundary[a] && !is_boundary[c]
                p[a] = (p[a] + p[c]) / 2
            end
            deg[a] += deg[c] - 2
            deg[best_b] -= 1
            deg[best_d] -= 1
            final_node[c] = a

            for j in 1:length(quads)
                is_dead[j] && continue
                replace!(quads[j], c => a)
                length(unique(quads[j])) < 4 && (is_dead[j] = true)
            end

        else  # :flip
            n0, n1, n2, n3, n4, n5 = best_flip_n
            steps = best_flip_steps
            if steps == 1
                # Shared edge rotates n0-n3 → n1-n4
                # quad A = [n0, n1, n4, n5],  quad B = [n1, n2, n3, n4]
                quads[best_flip_i1] = [n0, n1, n4, n5]
                quads[best_flip_i2] = [n1, n2, n3, n4]
                deg[n0] -= 1; deg[n3] -= 1
                deg[n1] += 1; deg[n4] += 1
            else  # steps == 2
                # Shared edge rotates n0-n3 → n2-n5
                # quad A = [n0, n1, n2, n5],  quad B = [n2, n3, n4, n5]
                quads[best_flip_i1] = [n0, n1, n2, n5]
                quads[best_flip_i2] = [n2, n3, n4, n5]
                deg[n0] -= 1; deg[n3] -= 1
                deg[n2] += 1; deg[n5] += 1
            end
        end

        live       = findall(.!is_dead)
        live_quads = [Index4(quads[j]...) for j in live]
        tmp        = DMesh(p, live_quads)
        adj        = [Int.(a) for a in node_adjacency(tmp)]
        D          = _node_depths(is_boundary, adj)
        nb_live    = element_face_neighbors(tmp)
    end

    live_quads = [Index4(quads[j]...) for j in live]

    # Path-compress final_node so chained merges (c→a, a→x) resolve to their root.
    for i in 1:np
        j = i
        while final_node[j] != j; j = final_node[j]; end
        final_node[i] = j
    end

    result = cleanup_mesh(DMesh(p, live_quads))

    # Build old→new node map.
    old_to_new = zeros(Int, np)
    for (new_id, old_id) in enumerate(result.ix)
        old_to_new[old_id] = new_id
    end
    # Merged/dead nodes: follow final_node to their surviving root.
    for i in 1:np
        old_to_new[i] == 0 && (old_to_new[i] = old_to_new[final_node[i]])
    end
    # Orphan nodes (in p but not referenced by any quad, e.g. triangle-only nodes):
    # append them to the output mesh so their IDs remain valid.
    np_out   = length(result.msh.p)
    orphan_p = eltype(p)[]
    for i in 1:np
        if old_to_new[i] == 0
            push!(orphan_p, p[i])
            old_to_new[i] = np_out + length(orphan_p)
        end
    end
    out_p = isempty(orphan_p) ? result.msh.p : vcat(result.msh.p, orphan_p)

    return DMesh(out_p, result.msh.t), old_to_new
end

################################################################################
# Public API
################################################################################

"""
    tri2quad(tmsh::DMesh, pfix; eps_rel=0.1, verbose=false) -> DMesh

Convert a triangular mesh to an all-quadrilateral mesh using a four-stage pipeline:

1. Topology-driven triangle matching via Blossom perfect matching
2. Boundary triangle collapse via union-find node merging
3. Greedy quad topology improvement (collapse + edge flip) to reduce degree irregularity (quads only; triangles held aside)
4. (Optional) Catmull-Clark refinement if triangles survive stages 2–3

`pfix` specifies the fixed corner node positions (same as passed to `distmesh2d`).
"""
function tri2quad(tmsh::DMesh, pfix=Point2d[]; eps_rel=0.1, verbose=false)
    # Stage 1: Match triangle pairs into quads
    topo = _QuadTopology(tmsh, pfix)
    q, t, tix = match_tri2quad(tmsh, topo)
    verbose && println("Stage 1: $(length(q)) quads, $(length(t)) unmatched triangles")

    # Stage 2: Collapse unmatched boundary triangles
    new_p, new_q, new_t = _collapse_boundary_triangles(tmsh, topo, q, t, tix)
    verbose && println("Stage 2: $(length(new_t)) triangle(s) survive after boundary collapse")

    # Stage 3: Greedy quad topology improvement (collapse + edge flip) on quads only;
    # surviving triangles are held aside.
    # Returns a node map so new_t can be remapped into the compacted mesh numbering.
    qmsh, node_map = _greedy_quad_improve(DMesh(new_p, new_q), pfix; verbose=verbose)
    verbose && println("Stage 3: $(length(qmsh.t)) quads after topology improvement")

    # Remap surviving triangles through the collapse node map.
    surviving_t = filter(rt -> length(unique(rt)) == 3,
                         [Index3(node_map[Int(tri[1])], node_map[Int(tri[2])], node_map[Int(tri[3])])
                          for tri in new_t])

    # Stage 4: Catmull-Clark refinement if triangles remain after collapse.
    if !isempty(surviving_t)
        @warn "$(length(surviving_t)) triangle(s) survived collapse; applying Catmull-Clark refinement"
        qmsh = _catmull_clark_refine(qmsh.p, qmsh.t, surviving_t)
    end

    return qmsh
end


function quad_project_and_smooth!(qmsh::DMesh, dfcn, hfcn, pfix=Point2d[];
                    plotting = false,          # Optional live plotting
                    maxiter = 10_000,          # When to terminate if no convergence
)
    deltat = 0.2
    h0 = minimum(norm.(all_edges(qmsh)))
    dptol = 1e-4 * h0
    deps = sqrt(eps()) * h0

    p = qmsh.p
    corner_ids = _corner_node_ids(p, pfix)
    boundary_ids = boundary_nodes(qmsh)
    project_ids = setdiff(boundary_ids, corner_ids)
    converged = false

    bars1 = all_edges(qmsh)
    bars2 = vcat([Index2(q[1], q[3]) for q in qmsh.t],
                 [Index2(q[2], q[4]) for q in qmsh.t])

    # Main loop
    for iter = 1:maxiter
        pold = copy(p)
        for (bars,Fsc) in ((bars1,0.8), (bars2,0.8))
            barvec = barvectors(p, bars)
            L = norm.(barvec)
            L0 = desiredlengths(p, bars, L, hfcn, Fsc)
            F = @. (L0^2 - L^2) / L
            p .+= deltat * total_node_forces(F, L, barvec, bars, length(p), corner_ids)
        end

        d = project_nodes!(p, dfcn, deps, project_ids)

        plotting && live_plot(qmsh)

        converged = maximum(norm.(p-pold); init=0.0) < dptol
        converged && break
    end
    
    converged || @warn "No convergence in maxiter=$maxiter iterations"

    return nothing
end

