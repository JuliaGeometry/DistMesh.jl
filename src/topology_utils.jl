################################################################################
### General mesh topology utilities
################################################################################

snap(x::T, scaling=1) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling*tol*round(x/scaling/tol) + zero(T)  # Adding zero to uniquify -0.0 and 0.0


"""
    cleanup_mesh(msh::DMesh) -> (msh::DMesh, ix::Vector{Int})

Remove duplicate nodes from the mesh `msh` and re-index the connectivity.

This function identifies nodes that are coincident (or within a very small tolerance relative to the mesh size) 
and merges them. This is useful after mesh generation or modification operations that might create 
overlapping vertices.

# Arguments
- `msh::DMesh`: The input mesh containing nodes `p` and connectivity `t`.

# Returns
A `NamedTuple` `(msh, ix)` where:
- `msh`: The cleaned `DMesh` with duplicate nodes removed.
- `ix`: An index vector mapping the **new** nodes to the **old** nodes (i.e., `new_p = old_p[ix]`).

# Example
```julia
clean_msh, = cleanup_mesh(dirty_msh)  # Ignoring the index output (ix)

```

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
    element_face_neighbors(msh::DMesh{D,T,E,N,I}) -> Matrix{Tuple{I, I}}

Compute element connectivities across faces.

Returns a matrix of size `(num_faces_per_element, num_elements)` where each entry 
is a tuple `(neighbor_element, neighbor_local_face)`.

If a face is on the boundary, the entry is `(0, 0)`.
"""
function element_face_neighbors(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    t = msh.t
    nt = length(t)
    
    fmap = facemap(E())
    nf = length(fmap)       
    nfv = length(fmap[1])   

    # Initialize a single matrix of tuples with (0, 0) for boundaries
    neighbors = fill((zero(I), zero(I)), nf, nt)
    
    # Key is SVector, Value is (element_idx, face_idx)
    dd = Dict{SVector{nfv, I}, Tuple{I, I}}()
    sizehint!(dd, nt * nf)

    for iel in 1:nt
        verts = t[iel]
        
        for jf in 1:nf
            key = sort(verts[fmap[jf]])
            
            if haskey(dd, key)
                nbel, nbface = pop!(dd, key)
                
                # Bi-directional link using tuples
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
    face_element_map(msh::DMesh{D,T,E,N,I}) -> Dict{SVector{nfv, I}, Vector{Tuple{I, I}}}

Compute the mapping from faces to all connected elements. 

Returns a dictionary where each key is a sorted `SVector` representing the face, 
and the value is a vector of tuples `(element_idx, local_face_idx)`.
"""
function face_element_map(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    t = msh.t
    nt = length(t)
    
    fmap = facemap(E())
    nf = length(fmap)       
    nfv = length(fmap[1])   

    # Dictionary maps sorted face -> Vector of (Element, Local Face)
    face_to_elements = Dict{SVector{nfv, I}, Vector{Tuple{I, I}}}()
    sizehint!(face_to_elements, div(nt * nf, 2))
    
    for iel in 1:nt
        verts = t[iel]
        
        for jf in 1:nf
            key = sort(verts[fmap[jf]])
            
            connected_elements = get!(face_to_elements, key, Tuple{I, I}[])
            push!(connected_elements, (iel, jf))
        end
    end
    
    return face_to_elements
end

"""
    filter_elements_by_degree(face_map, condition)

Core function to extract unique element indices from a face map where the number 
of elements sharing a face satisfies the given `condition` function.
"""
function filter_elements_by_degree(face_map, condition)
    problem_elements = [first(item) for v in values(face_map) if condition(length(v)) for item in v]
    return unique!(problem_elements)
end

find_nonmanifold_elements(face_map) = filter_elements_by_degree(face_map, x -> x > 2)
find_boundary_elements(face_map)    = filter_elements_by_degree(face_map, x -> x == 1)

"""
    is_manifold_mesh(m::DMesh) -> Bool

Check if the mesh is manifold.

A mesh is manifold if every face (edge in 2D, triangle in 3D) is shared by exactly 2 elements.
This is equivalent to checking that the mesh has no non-manifold features (edges/faces touching more than 2 elements).

# Arguments
- `m`: The mesh to check.

# Returns
- `true` if the mesh is manifold, `false` otherwise.
"""
function is_manifold_mesh(m::DMesh)
    e2f = face_element_map(m)
    return all(length.(values(e2f)) .== 2)
end

"""
    foreach_face(f::Function, msh::DMesh)

Iterate over all local faces of all elements in the mesh, applying the function `f`.

This is an internal higher-order helper function designed to traverse the mesh 
topology efficiently. It handles the boilerplate of looking up element-to-element 
connectivity and face mappings.

The provided function `f` must accept four arguments:
1. `iel`: The index of the current element.
2. `jf`: The local index of the face within the current element.
3. `jel`: The index of the neighboring element sharing this face (0 if it is a boundary face).
4. `map`: The local node mapping for the faces of this element type.

# Example Usage
```julia
foreach_face(msh) do iel, jf, jel, map
    if jel == 0
        println("Found a boundary face on element \$iel")
    end
end
"""
function foreach_face(f::Function, msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    nb = element_face_neighbors(msh)
    nt = length(msh.t)
    map = facemap(E())
    nf = length(map)

    for iel in 1:nt
        for jf in 1:nf
            jel = nb[jf, iel][1]
            # Call the user-provided function with the current state
            f(iel, jf, jel, map)
        end
    end
end

"""
    boundary_faces(msh::DMesh{D,T,N,I}) -> Vector{SVector{L, I}}

Identify the boundary faces of the mesh.

Returns a list of all mesh faces that are not shared by two elements. 
- For a 2D triangular mesh, these are the boundary edges (line segments).
- For a 3D tetrahedral mesh, these are the boundary faces (triangles).

# Arguments
- `msh`: The mesh object.

# Returns
- A `Vector` of `SVector`s, where each `SVector` contains the node indices of a boundary face.
"""
function boundary_faces(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    map = facemap(E())
    nfv = length(map[1])
    nt = length(msh.t)
    
    bnd = SVector{nfv,I}[]
    sizehint!(bnd, floor(Int, sqrt(nt) * 4)) 

    foreach_face(msh) do iel, jf, jel, map
        if jel == 0
            push!(bnd, msh.t[iel][map[jf]])
        end
    end

    return bnd
end

"""
    all_faces(msh::DMesh{D,T,N,I}) -> Tuple{Vector{SVector{L, I}}, Vector{Int}}

Identify all unique faces in the mesh and locate the boundary faces.

Returns a complete list of every face in the mesh exactly once. For interior faces 
shared by two adjacent elements, only a single instance is recorded. Additionally, 
it returns the indices of the faces that lie on the external boundary (i.e., faces 
not shared by another element).

- For a 2D triangular mesh, these are all unique edges in the mesh.
- For a 3D tetrahedral mesh, these are all unique triangular faces.

# Arguments
- `msh`: The mesh object.

# Returns
- A `Tuple` containing two vectors:
    1. `faces`: A `Vector` of `SVector`s, where each `SVector` contains the node indices of a unique face.
    2. `boundary_idx`: A `Vector{Int}` containing the corresponding indices of the boundary faces within the `faces` array.
"""
function all_faces(msh::DMesh{D,T,E,N,I}) where {D,T,E,N,I}
    nfv = length(facemap(E())[1])
    
    faces = SVector{nfv,I}[]
    boundary_idx = Int[]

    foreach_face(msh) do iel, jf, jel, map
        if jel == 0 || jel > iel
            push!(faces, msh.t[iel][map[jf]])
        end
        if jel == 0
            push!(boundary_idx, length(faces))
        end
    end

    return faces, boundary_idx
end

"""
    boundary_nodes(msh::DMesh) -> Vector{I}

Identify the boundary nodes (vertices) of the mesh.

Returns a list of all unique node indices that lie on the external boundary 
of the mesh. These are the nodes that make up the boundary faces (in 3D) 
or boundary edges (in 2D).

# Arguments
- `msh`: The mesh object.

# Returns
- A `Vector` of integers (of type `I`) containing the unique indices of all boundary nodes.
"""
boundary_nodes(msh::DMesh) = unique(Iterators.flatten(boundary_faces(msh)))

"""
    all_edges(msh::DMesh) -> Vector{SVector{2, I}}

Identify all unique edges in the mesh.

Extracts all edges from every element, normalizes their orientation (smallest node index first), 
and returns a list of unique edges. Works for both 2D and 3D meshes.

# Arguments
- `msh`: The mesh object.

# Returns
- A `Vector` of 2-element `SVector`s representing the unique edges in the mesh.
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

"""
    node_degrees(msh::DMesh{2}) -> Vector{Int}

Compute the degree of each node (vertex) in a 2D mesh.

The degree is calculated as the number of unique edges connected to a given node. 

# Arguments
- `msh`: A 2D mesh object (`DMesh{2}`).

# Returns
- A `Vector{Int}` of the same length as the number of nodes in the mesh, 
  where the `i`-th entry contains the degree of the `i`-th node.
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

Compute the adjacency list of each node (vertex) in the mesh.

Returns a vector of length equal to the number of nodes, where the `i`-th entry
is a sorted list of the indices of all nodes connected to node `i` by an edge.

The adjacency is derived directly from the element edge map, avoiding the
intermediate sorted edge array that `all_edges` produces.

# Arguments
- `msh`: The mesh object.

# Returns
- A `Vector{Vector{I}}` where `I` is the mesh index type, of length `length(msh.p)`.
  Each inner vector contains the sorted indices of the neighboring nodes.
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
