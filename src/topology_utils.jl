################################################################################
### Internal helpers
################################################################################

snap(x::T, scaling=1) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling*tol*round(x/scaling/tol) + zero(T)  # Adding zero to uniquify -0.0 and 0.0

function filter_elements_by_degree(face_map, condition)
    elems = [first(item) for v in values(face_map) if condition(length(v)) for item in v]
    return unique!(elems)
end

"""
    foreach_face(f::Function, msh::DMesh)

Iterate over every local face of every element, calling `f(iel, jf, jel, fmap)`.

- `iel`: element index
- `jf`: local face index within the element
- `jel`: index of the neighboring element (0 for boundary faces)
- `fmap`: face-to-local-node map for the element type
"""
function foreach_face(f::Function, msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    nb = element_face_neighbors(msh)
    nt = length(msh.t)
    fmap = facemap(E())
    nf = length(fmap)

    for iel in 1:nt
        for jf in 1:nf
            jel = nb[jf, iel][1]
            f(iel, jf, jel, fmap)
        end
    end
end

################################################################################
### Mesh utilities
################################################################################

"""
    cleanup_mesh(msh::DMesh) -> (msh::DMesh, ix::Vector{Int})

Remove duplicate nodes and re-index the connectivity.

Nodes that are coincident (within a small tolerance relative to the mesh size)
are merged. Returns a named tuple `(msh, ix)` where `ix` maps new node indices
to old ones (`new_p = old_p[ix]`).
"""
function cleanup_mesh(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    p, t = msh

    scaling = maximum(norm.(p))
    scaling == 0.0 && (scaling = 1.0)
    pp = [snap.(p1, scaling) for p1 in p]
    ppp = unique(pp)
    ix = I.(indexin(ppp, pp))
    jx = I.(indexin(pp, ppp))

    tt = jx[reinterpret(I, t)]
    pix = unique(tt)
    jx1 = I.(indexin(tt, pix))

    new_p = p[ix[pix]]
    new_t = collect(reinterpret(SVector{N, I}, jx1))
    ix_final = ix[pix]

    return (msh=DMesh(new_p, new_t), ix=ix_final)
end

"""
    is_manifold_mesh(msh::DMesh) -> Bool

Return `true` if every face is shared by exactly two elements, `false` otherwise.
"""
function is_manifold_mesh(msh::DMesh)
    e2f = face_element_map(msh)
    return all(length.(values(e2f)) .== 2)
end

################################################################################
### Element utilities
################################################################################

"""
    element_face_neighbors(msh::DMesh) -> Matrix{Tuple{I, I}}

Compute element-to-element connectivity across faces.

Returns a `(num_faces_per_element, num_elements)` matrix where each entry is a
tuple `(neighbor_element, neighbor_local_face)`. Boundary faces have entry `(0, 0)`.
"""
function element_face_neighbors(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    t = msh.t
    nt = length(t)

    fmap = facemap(E())
    nf = length(fmap)
    nfv = length(fmap[1])

    neighbors = fill((zero(I), zero(I)), nf, nt)
    dd = Dict{SVector{nfv, I}, Tuple{I, I}}()
    sizehint!(dd, nt * nf)

    for iel in 1:nt
        verts = t[iel]
        for jf in 1:nf
            key = sort(verts[fmap[jf]])
            if haskey(dd, key)
                nbel, nbface = pop!(dd, key)
                neighbors[jf, iel] = (nbel, nbface)
                neighbors[nbface, nbel] = (iel, jf)
            else
                dd[key] = (iel, jf)
            end
        end
    end

    return neighbors
end

"""
    face_element_map(msh::DMesh) -> Dict{SVector, Vector{Tuple{I, I}}}

Map each face to the elements that contain it.

Returns a dictionary keyed by sorted face node indices; each value is a vector of
`(element_idx, local_face_idx)` tuples. Interior faces have two entries; boundary
faces have one.
"""
function face_element_map(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    t = msh.t
    nt = length(t)

    fmap = facemap(E())
    nf = length(fmap)
    nfv = length(fmap[1])

    face_to_elements = Dict{SVector{nfv, I}, Vector{Tuple{I, I}}}()
    sizehint!(face_to_elements, div(nt * nf, 2))

    for iel in 1:nt
        verts = t[iel]
        for jf in 1:nf
            key = sort(verts[fmap[jf]])
            push!(get!(face_to_elements, key, Tuple{I, I}[]), (iel, jf))
        end
    end

    return face_to_elements
end

"""
    find_boundary_elements(msh::DMesh) -> Vector{I}

Return the indices of elements that have at least one boundary face.
"""
find_boundary_elements(msh::DMesh) =
    filter_elements_by_degree(face_element_map(msh), x -> x == 1)

"""
    find_nonmanifold_elements(msh::DMesh) -> Vector{I}

Return the indices of elements that share a face with more than one other element.
"""
find_nonmanifold_elements(msh::DMesh) =
    filter_elements_by_degree(face_element_map(msh), x -> x > 2)

################################################################################
### Face and edge utilities
################################################################################

"""
    all_faces(msh::DMesh) -> (Vector{SVector{L,I}}, Vector{Int})

Return all unique faces in the mesh together with the indices of the boundary faces.

Each interior face is listed once; boundary faces are listed once and their positions
in the output vector are recorded in `boundary_idx`.
For 2D triangular meshes these are edges; for 3D tetrahedral meshes, triangles.
"""
function all_faces(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    nfv = length(facemap(E())[1])

    faces = SVector{nfv,I}[]
    boundary_idx = Int[]

    foreach_face(msh) do iel, jf, jel, fmap
        if jel == 0 || jel > iel
            push!(faces, msh.t[iel][fmap[jf]])
        end
        if jel == 0
            push!(boundary_idx, length(faces))
        end
    end

    return faces, boundary_idx
end

"""
    boundary_faces(msh::DMesh) -> Vector{SVector{L,I}}

Return all faces on the boundary of the mesh (faces not shared by two elements).

For 2D triangular meshes these are boundary edges; for 3D tetrahedral meshes,
boundary triangles.
"""
function boundary_faces(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    fmap = facemap(E())
    nfv = length(fmap[1])
    nt = length(msh.t)

    bnd = SVector{nfv,I}[]
    sizehint!(bnd, floor(Int, sqrt(nt) * 4))

    foreach_face(msh) do iel, jf, jel, fmap
        if jel == 0
            push!(bnd, msh.t[iel][fmap[jf]])
        end
    end

    return bnd
end

"""
    all_edges(msh::DMesh) -> Vector{SVector{2,I}}

Return all unique edges in the mesh, each oriented with the smaller node index first.
"""
function all_edges(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    emap = edgemap(E())

    total_edges = length(msh.t) * length(emap)
    edges = Vector{SVector{2, I}}(undef, total_edges)

    idx = 1
    for el in msh.t
        for e_local in emap
            n1, n2 = el[e_local]
            edges[idx] = SVector(min(n1, n2), max(n1, n2))
            idx += 1
        end
    end

    return unique!(sort!(edges))
end

################################################################################
### Node utilities
################################################################################

"""
    boundary_nodes(msh::DMesh) -> Vector{I}

Return the indices of all nodes that lie on the boundary of the mesh.
"""
boundary_nodes(msh::DMesh) = unique(Iterators.flatten(boundary_faces(msh)))

"""
    node_degrees(msh::DMesh) -> Vector{Int}

Return the degree (number of incident edges) of each node.
"""
function node_degrees(msh::DMesh)
    edges = all_edges(msh)
    deg = zeros(Int, length(msh.p))

    for e in edges
        for i in e
            deg[i] += 1
        end
    end

    return deg
end

"""
    node_adjacency(msh::DMesh) -> Vector{Vector{I}}

Return the adjacency list of each node.

The `i`-th entry is a sorted vector of indices of all nodes connected to node `i`
by an edge.
"""
function node_adjacency(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    emap = edgemap(E())
    adj = [I[] for _ in 1:length(msh.p)]

    for el in msh.t
        for e_local in emap
            n1, n2 = el[e_local]
            push!(adj[n1], n2)
            push!(adj[n2], n1)
        end
    end

    for a in adj
        sort!(unique!(a))
    end

    return adj
end

"""
    node_element_map(msh::DMesh) -> Vector{Vector{I}}

Return the list of elements containing each node.

The `i`-th entry is a sorted vector of indices of all elements that include node `i`.
"""
function node_element_map(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    np = length(msh.p)
    nt = length(msh.t)
    nemap = [I[] for _ in 1:np]

    for iel in 1:nt
        for n in msh.t[iel]
            push!(nemap[n], I(iel))
        end
    end

    for a in nemap
        sort!(a)
    end

    return nemap
end

"""
    segcollect(bedges::Vector{SVector{2,I}}) -> Vector{Vector{I}}

Stitch unordered boundary edges into ordered, closed boundary loops.

Each returned segment is a vector of node indices forming a closed loop
(the first node is NOT repeated at the end). The orientation follows the
natural ordering of the boundary edges as returned by `boundary_faces`.

This is useful for computing interior angles at boundary nodes or for
identifying corners in the boundary geometry.
"""
function segcollect(bedges::Vector{SVector{2,I}}) where {I <: Integer}
    # Build adjacency map: node -> [neighbor1, neighbor2] (boundary edges only)
    adj = Dict{I, Vector{I}}()
    for e in bedges
        push!(get!(adj, e[1], I[]), e[2])
        push!(get!(adj, e[2], I[]), e[1])
    end

    visited_edges = Set{Tuple{I,I}}()
    segments = Vector{Vector{I}}()

    for e in bedges
        key = (min(e[1], e[2]), max(e[1], e[2]))
        key in visited_edges && continue

        # Start a new segment from e[1] → e[2]
        seg = I[e[1]]
        push!(visited_edges, key)
        curr = e[2]
        prev = e[1]

        while curr != seg[1]
            push!(seg, curr)
            nbrs = adj[curr]
            # Pick the neighbor that isn't where we came from
            nxt = nbrs[1] == prev ? nbrs[2] : nbrs[1]
            key2 = (min(curr, nxt), max(curr, nxt))
            push!(visited_edges, key2)
            prev = curr
            curr = nxt
        end

        push!(segments, seg)
    end

    return segments
end
