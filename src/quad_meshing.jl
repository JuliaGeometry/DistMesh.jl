################################################################################
# quad_meshing.jl  –  Triangle-to-quad conversion for DistMesh.jl
#
# Pipeline:
#   Stage 1: Topology-driven tri-to-quad matching (Blossom perfect matching)
#   Stage 2: Boundary triangle collapse (union-find node merging)
#   Stage 3: Optional Catmull-Clark refinement (if triangles survive stage 2)
#   Stage 4: Greedy quad collapse (depth-weighted diagonal collapse)
#   Stage 5: Dart untangling (analytical Jacobian fix for inverted quads)
################################################################################

################################################################################
# Internal topology cache
################################################################################

struct _QuadTopology
    is_corner      :: Vector{Bool}
    is_boundary    :: Vector{Bool}
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
function _element_depths(msh::DMesh)
    nb = element_face_neighbors(msh)
    nt = length(msh.t)
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
            min_d[i]   = 2
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

    adj = [Int.(a) for a in node_adjacency(msh)]
    deg = node_degrees(msh)
    ideal_d, min_d = _target_degrees(msh, Vector{Bool}(is_boundary), Vector{Bool}(is_corner))

    return _QuadTopology(
        Vector{Bool}(is_corner),
        Vector{Bool}(is_boundary),
        _element_depths(msh),
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

        gain_u = 2 * (d_curr[u] - d_ideal[u]) - 1
        gain_v = 2 * (d_curr[v] - d_ideal[v]) - 1
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
function match_tri2quad(msh::DMesh, pfix)
    topo = _QuadTopology(msh, pfix)
    nb   = element_face_neighbors(msh)
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

function _collapse_boundary_triangles(msh::DMesh, q::Vector{Index4},
                                       t_unmatched::Vector{<:SVector{3}},
                                       t_global_idx::Vector{Int},
                                       is_corner::Vector{Bool})
    np  = length(msh.p)
    nb  = element_face_neighbors(msh)
    emap = edgemap(Simplex{2}())

    parent      = collect(1:np)
    corner_root = zeros(Int, np)
    for i in 1:np
        is_corner[i] && (corner_root[i] = i)
    end

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

    # Merge the two nodes on the first boundary edge of each unmatched triangle
    for gt in t_global_idx
        e = findfirst(ie -> Int(nb[ie, gt][1]) == 0, 1:3)
        e !== nothing && union!(Int(msh.t[gt][emap[e][1]]), Int(msh.t[gt][emap[e][2]]))
    end

    # Path-compress
    for i in 1:np
        parent[i] = find_root(i)
    end

    # Build old→new ID map and new coordinates
    unique_roots = unique(parent)
    old_to_new   = zeros(Int, np)
    new_p        = Vector{Point2d}(undef, length(unique_roots))

    for (new_id, root) in enumerate(unique_roots)
        members = findall(==(root), parent)
        for old_id in members
            old_to_new[old_id] = new_id
        end
        if corner_root[root] != 0
            new_p[new_id] = msh.p[corner_root[root]]
        else
            new_p[new_id] = sum(msh.p[m] for m in members) / length(members)
        end
    end

    # Remap quads; drop degenerate (repeated-node) quads
    new_q = Index4[]
    for quad in q
        rq = Index4(old_to_new[quad[1]], old_to_new[quad[2]],
                    old_to_new[quad[3]], old_to_new[quad[4]])
        length(unique(rq)) == 4 && push!(new_q, rq)
    end

    # Remap surviving triangles; drop degenerate
    new_t = Index3[]
    for tri in t_unmatched
        rt = Index3(old_to_new[tri[1]], old_to_new[tri[2]], old_to_new[tri[3]])
        length(unique(rt)) == 3 && push!(new_t, rt)
    end

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
# Stage 4: Greedy quad collapse
################################################################################

function _delta_E(a, b, c, d, deg, ideal_deg, D)
    W(n) = 100 - D[n]
    E(n, dnew) = W(n) * (dnew - ideal_deg[n])^2
    E_old = E(a, deg[a]) + E(b, deg[b]) + E(c, deg[c]) + E(d, deg[d])
    E_new = E(a, deg[a] + deg[c] - 2) + E(b, deg[b] - 1) + E(d, deg[d] - 1)
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

function _greedy_quad_collapse(msh::DMesh, pfix; verbose=false)
    topo = _QuadTopology(msh, pfix)

    p         = copy(msh.p)
    quads     = [[Int(n) for n in q] for q in msh.t]  # mutable Vector{Vector{Int}}
    is_dead   = falses(length(quads))
    deg       = copy(topo.degree)
    ideal_deg = topo.ideal_degree
    is_boundary = topo.is_boundary
    is_corner   = topo.is_corner
    adj = [Int.(a) for a in node_adjacency(msh)]
    D   = _node_depths(is_boundary, adj)

    for pass in 1:100
        best_dE, best_a, best_c, best_b, best_d = 0.0, 0, 0, 0, 0

        for i in 1:length(quads)
            is_dead[i] && continue
            ns = quads[i]
            for (ai, bi, ci, di) in ((1,2,3,4), (2,3,4,1))
                a, b, c, d = ns[ai], ns[bi], ns[ci], ns[di]
                _collapse_safe(a, c, is_boundary, is_corner, adj) || continue
                dE = _delta_E(a, b, c, d, deg, ideal_deg, D)
                if dE < best_dE
                    best_dE = dE
                    best_a = is_boundary[a] ? a : (is_boundary[c] ? c : a)
                    best_c = (best_a == a) ? c : a
                    best_b, best_d = b, d
                end
            end
        end

        best_dE >= -0.01 && (verbose && println("  Quad collapse converged in $pass passes."); break)

        a, c = best_a, best_c
        if !is_boundary[a] && !is_boundary[c]
            p[a] = (p[a] + p[c]) / 2
        end
        deg[a] += deg[c] - 2
        deg[best_b] -= 1
        deg[best_d] -= 1

        for j in 1:length(quads)
            is_dead[j] && continue
            replace!(quads[j], c => a)
            length(unique(quads[j])) < 4 && (is_dead[j] = true)
        end

        live = findall(.!is_dead)
        live_quads = [Index4(quads[j]...) for j in live]
        tmp = DMesh(p, live_quads)
        adj = [Int.(a) for a in node_adjacency(tmp)]
        D   = _node_depths(is_boundary, adj)
    end

    live = findall(.!is_dead)
    live_quads = [Index4(quads[j]...) for j in live]
    return cleanup_mesh(DMesh(p, live_quads)).msh
end

################################################################################
# Stage 5: Dart untangling
################################################################################

function _untangle_darts!(p::Vector{Point2d}, q::Vector{Index4},
                           is_boundary::Vector{Bool}, is_corner::Vector{Bool};
                           eps_rel=1e-1)
    fixes = 0

    calc_J(a, b, d) = (p[b][1]-p[a][1])*(p[d][2]-p[a][2]) -
                      (p[b][2]-p[a][2])*(p[d][1]-p[a][1])
    dist_sq(n1, n2) = (p[n2][1]-p[n1][1])^2 + (p[n2][2]-p[n1][2])^2

    for quad in q
        n1, n2, n3, n4 = Int(quad[1]), Int(quad[2]), Int(quad[3]), Int(quad[4])
        J1 = calc_J(n1, n2, n4)
        J2 = calc_J(n2, n3, n1)
        J3 = calc_J(n3, n4, n2)
        J4 = calc_J(n4, n1, n3)

        min_J, min_idx = findmin(SA[J1, J2, J3, J4])
        a, b, _, d = if min_idx == 1; (n1, n2, n3, n4)
                     elseif min_idx == 2; (n2, n3, n4, n1)
                     elseif min_idx == 3; (n3, n4, n1, n2)
                     else                 (n4, n1, n2, n3)
                     end

        eps_abs = eps_rel * 0.5 * (dist_sq(a, b) + dist_sq(a, d))
        min_J > eps_abs - 1e-8 && continue
        fixes += 1

        if !is_boundary[a]
            gx = p[b][2] - p[d][2]
            gy = p[d][1] - p[b][1]
            norm_sq = gx^2 + gy^2
            if norm_sq > 1e-12
                step = (eps_abs - min_J) / norm_sq
                p[a] = p[a] + Point2d(step * gx, step * gy)
            end
        elseif is_corner[a]
            # Corner locked: move an interior neighbor instead
            move_node = !is_boundary[b] ? b : (!is_boundary[d] ? d : 0)
            if move_node != 0
                if move_node == b
                    gx = p[d][2] - p[a][2];  gy = p[a][1] - p[d][1]
                else
                    gx = p[a][2] - p[b][2];  gy = p[b][1] - p[a][1]
                end
                norm_sq = gx^2 + gy^2
                if norm_sq > 1e-12
                    step = (eps_abs - min_J) / norm_sq
                    p[move_node] = p[move_node] + Point2d(step * gx, step * gy)
                end
            end
        end
        # Non-corner boundary reflex: skip silently (should be rare after smoothing)
    end

    return fixes
end

function _untangle_darts_sweep!(p::Vector{Point2d}, q::Vector{Index4},
                                 is_boundary::Vector{Bool}, is_corner::Vector{Bool};
                                 eps_rel=1e-1, verbose=false)
    for sweep in 1:10
        fixes = _untangle_darts!(p, q, is_boundary, is_corner; eps_rel=eps_rel)
        if fixes == 0
            verbose && println("  Untangled in $sweep sweep(s).")
            return
        end
    end
    verbose && @warn "Untangler reached max sweeps; mesh may still contain inverted elements."
end

################################################################################
# Public API
################################################################################

"""
    tri2quad(tmsh::DMesh, pfix; eps_rel=0.1, verbose=false) -> DMesh

Convert a triangular mesh to an all-quadrilateral mesh using a four-stage pipeline:

1. Topology-driven triangle matching via Blossom perfect matching
2. Boundary triangle collapse via union-find node merging
3. (Optional) Catmull-Clark refinement if triangles survive stage 2
4. Greedy quad collapse to reduce degree irregularity
5. Analytical dart untangling to fix inverted quads

`pfix` specifies the fixed corner node positions (same as passed to `distmesh2d`).
"""
function tri2quad(tmsh::DMesh, pfix; eps_rel=0.1, verbose=false)
    # Stage 1: Match triangle pairs into quads
    q, t, tix = match_tri2quad(tmsh, pfix)
    verbose && println("Stage 1: $(length(q)) quads, $(length(t)) unmatched triangles")

    # Stage 2: Collapse unmatched boundary triangles
    topo0 = _QuadTopology(tmsh, pfix)
    new_p, new_q, new_t = _collapse_boundary_triangles(tmsh, q, t, tix, topo0.is_corner)
    verbose && println("Stage 2: $(length(new_t)) triangle(s) survive after boundary collapse")
#    return DMesh(new_p, new_t), DMesh(new_p, new_q)

    # Stage 3: Catmull-Clark refinement if triangles remain
    if !isempty(new_t)
        verbose && @warn "$(length(new_t)) triangle(s) survived collapse; applying Catmull-Clark refinement"
        qmsh = _catmull_clark_refine(new_p, new_q, new_t)
    else
        qmsh = DMesh(new_p, new_q)
    end

    # Stage 4: Greedy quad collapse
    qmsh = _greedy_quad_collapse(qmsh, pfix; verbose=verbose)
    verbose && println("Stage 4: $(length(qmsh.t)) quads after collapse")

    # Stage 5: Untangle inverted quads
    topo_q = _QuadTopology(qmsh, pfix)
    p_mut  = copy(qmsh.p)
    _untangle_darts_sweep!(p_mut, qmsh.t, topo_q.is_boundary, topo_q.is_corner;
                            eps_rel=eps_rel, verbose=verbose)
    return DMesh(p_mut, qmsh.t)
end

function quad_project_and_smooth!(qmsh::DMesh, fd, fh, pfix)
    @warn "quad_project_and_smooth! is not yet implemented - doing nothing"
end
