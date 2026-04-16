################################################################################
### Delaunator wrappers
################################################################################

extract_elements(tri::Delaunator.Triangulation{I}) where {I} =
    reinterpret(SVector{3, I}, triangles(tri))

function fix_winding!(tri::Delaunator.Triangulation{I}) where {I}
    tris = triangles(tri)
    mat = reinterpret(reshape, I, tris)
    
    # Swap rows 2 and 3 (in-place)
    @inbounds for c in axes(mat, 2)
        mat[2, c], mat[3, c] = mat[3, c], mat[2, c]
    end
    return nothing
end

function delaunay(p)
    t = triangulate(p)
    fix_winding!(t)
    return extract_elements(t)
end

################################################################################
### Size function utilities
################################################################################

"""
    huniform(p)

Returns `1.0`. Default sizing function for uniform meshes.
"""
huniform(p) = 1

################################################################################
### Element Mapping / Generators
################################################################################

"""
    element_map(f, msh::DMesh{D, T, G})

Return a generator that lazily applies `f(G(), nodes)` to the nodes of each 
element in the mesh. `nodes` is an `SVector` of the physical coordinates.

Example:
`total_vol = sum(element_map(element_volume, msh))`
"""
function element_map(f, msh::DMesh{D, T, G}) where {D, T, G}
    return (f(G(), msh.p[el]) for el in msh.t)
end

################################################################################
### Element Volumes
################################################################################

element_volume(::ElementGeometry, el) = error("Not implemented for this geometry")

"""
    element_volume(::Simplex{1}, el)

Compute the length of a 1D line segment.
"""
element_volume(::Simplex{1}, el) = norm(el[2] - el[1])

"""
    element_volume(::Block{1}, el)

Compute the length of a 1D line segment.
"""
element_volume(::Block{1}, el) = norm(el[2] - el[1])

"""
    element_volume(::Simplex{2}, el)

Compute the area of a 2D triangle.
"""
function element_volume(::Simplex{2}, el)
    p1, p2, p3 = el
    p12 = p2 - p1
    p13 = p3 - p1
    return (p12[1] * p13[2] - p12[2] * p13[1]) / 2
end

"""
    element_volume(::Simplex{3}, el)

Compute the volume of a 3D tetrahedron.
"""
function element_volume(::Simplex{3}, el)
    p1, p2, p3, p4 = el
    return dot(p2 - p1, cross(p3 - p1, p4 - p1)) / 6
end

"""
    element_volume(::Block{2}, el)

Compute the area of a 2D quadrilateral.
"""
function element_volume(::Block{2}, el)
    p1,p2,p3,p4 = el
    AC = p3 - p1
    BD = p4 - p2
    return 0.5*(AC[1]*BD[2] - AC[2]*BD[1])
end

"""
    element_volume(::Block{3}, el)

Compute the volume of a 3D hexahedron.
"""
function element_volume(::Block{3}, el)
    p1,p2,p3,p4,p5,p6,p7,p8 = el
    # Decomposes into 5 tetrahedra.
    tet(a,b,c,d) = dot(b-a, cross(c-a, d-a)) / 6
    return tet(p1,p2,p4,p5) + tet(p2,p3,p4,p7) + tet(p2,p5,p6,p7) +
           tet(p4,p5,p7,p8) + tet(p2,p4,p5,p7)
end

################################################################################
### Element Qualities - Simplex Elements
################################################################################

# Metric 1: Radius Ratio
element_quality_radius_ratio(::Simplex{1}, el) = 1.0

function element_quality_radius_ratio(::Simplex{2}, el)
    p1, p2, p3 = el
    a = norm(p2 - p1)
    b = norm(p3 - p2)
    c = norm(p1 - p3)
    
    s = (a + b + c) / 2
    denom = a * b * c
    
    return denom ≈ 0 ? 0.0 : 8 * (s - a) * (s - b) * (s - c) / denom
end

function element_quality_radius_ratio(::Simplex{3}, el)
    p1, p2, p3, p4 = el
    d12 = p2 - p1;  d13 = p3 - p1;  d14 = p4 - p1
    d23 = p3 - p2;  d24 = p4 - p2;  d34 = p4 - p3

    v  = element_volume(Simplex{3}(), el)
    s1 = norm(d12 × d13) / 2
    s2 = norm(d12 × d14) / 2
    s3 = norm(d13 × d14) / 2
    s4 = norm(d23 × d24) / 2

    p_1 = norm(d12) * norm(d34)
    p_2 = norm(d23) * norm(d14)
    p_3 = norm(d13) * norm(d24)

    denom = (s1 + s2 + s3 + s4) * sqrt((p_1 + p_2 + p_3) * (p_1 + p_2 - p_3) *
                                        (p_1 + p_3 - p_2) * (p_2 + p_3 - p_1))
    return denom ≈ 0 ? 0.0 : 216 * v^2 / denom
end

# Metric 2: Mean Ratio
element_quality_mean_ratio(::Simplex{1}, el) = 1.0

function element_quality_mean_ratio(::Simplex{2}, el)
    p1, p2, p3 = el
    area = element_volume(Simplex{2}(), el)
    l_sq = sum(abs2, p2 - p1) + sum(abs2, p3 - p2) + sum(abs2, p1 - p3)
    
    return l_sq ≈ 0 ? 0.0 : (4 * sqrt(3) * area) / l_sq
end

function element_quality_mean_ratio(::Simplex{3}, el)
    p1, p2, p3, p4 = el
    v    = element_volume(Simplex{3}(), el)
    l_sq = sum(abs2, p2 - p1) + sum(abs2, p3 - p1) + sum(abs2, p4 - p1) +
           sum(abs2, p3 - p2) + sum(abs2, p4 - p2) + sum(abs2, p4 - p3)

    return l_sq ≈ 0 ? 0.0 : 216 * v / sqrt(3) / l_sq^(3/2)
end

################################################################################
### Element Qualities - Block Elements
################################################################################

# Mean ratio: normalized to [0,1], 1 = perfect square/cube. Analogous to
# element_quality_mean_ratio for simplices (Knupp 2000).
function element_quality_mean_ratio(::Block{2}, el)
    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4)  # edge vectors (cyclic)
    cross2d(u, v) = u[1]*v[2] - u[2]*v[1]

    l2 = SVector(sum(abs2, e[1]), sum(abs2, e[2]), sum(abs2, e[3]), sum(abs2, e[4]))
    A  = SVector(cross2d(e[1],-e[4]), cross2d(e[2],-e[1]),
                 cross2d(e[3],-e[2]), cross2d(e[4],-e[3]))
    Q  = SVector((l2[1]+l2[3]), (l2[2]+l2[4]),
                 (l2[3]+l2[1]), (l2[4]+l2[2])) ./ (2 .* A)

    return 4 / sum(Q)
end

function element_quality_mean_ratio(::Block{3}, el)
    p1,p2,p3,p4,p5,p6,p7,p8 = el
    # At each corner, form the 3x3 Jacobian W from the 3 incident edge vectors.
    # Corner ordering matches the hex node layout (bottom 1-2-3-4, top 5-6-7-8).
    corners = (
        (p2-p1, p4-p1, p5-p1),  # corner 1
        (p3-p2, p1-p2, p6-p2),  # corner 2
        (p4-p3, p2-p3, p7-p3),  # corner 3
        (p1-p4, p3-p4, p8-p4),  # corner 4
        (p8-p5, p6-p5, p1-p5),  # corner 5
        (p5-p6, p7-p6, p2-p6),  # corner 6
        (p6-p7, p8-p7, p3-p7),  # corner 7
        (p7-p8, p5-p8, p4-p8),  # corner 8
    )
    function corner_quality(e1, e2, e3)
        detW = dot(e1, cross(e2, e3))
        detW <= 0 && return 0.0
        frob2 = sum(abs2, e1) + sum(abs2, e2) + sum(abs2, e3)
        return 3 * cbrt(detW^2) / frob2
    end
    return minimum(corner_quality(c...) for c in corners)
end

# Condition number of the corner Jacobian matrix (Knupp 2000).
function element_quality_condition_number(::Block{2}, el)
    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4)  # edge vectors (cyclic)
    cross2d(u, v) = u[1]*v[2] - u[2]*v[1]

    l2   = SVector(sum(abs2, e[1]), sum(abs2, e[2]), sum(abs2, e[3]), sum(abs2, e[4]))
    sins = SVector(cross2d(e[1],-e[4]), cross2d(e[2],-e[1]),
                   cross2d(e[3],-e[2]), cross2d(e[4],-e[3])) ./
           SVector(sqrt(l2[4]*l2[1]), sqrt(l2[1]*l2[2]),
                   sqrt(l2[2]*l2[3]), sqrt(l2[3]*l2[4]))
    any(<=(0), sins) && return 0.0

    k = SVector(l2[4]+l2[1], l2[1]+l2[2], l2[2]+l2[3], l2[3]+l2[4]) ./
        (SVector(sqrt(l2[4]*l2[1]), sqrt(l2[1]*l2[2]),
                 sqrt(l2[2]*l2[3]), sqrt(l2[3]*l2[4])) .* sins)

    return 4 / sqrt(sum(abs2, k))
end

function element_quality_min_scaled_jacobian(::Block{2}, el)
    # Minimum scaled corner Jacobian in [-1,1]; <= 0 means concave or self-intersecting
    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4)  # edge vectors (cyclic)
    cross2d(u, v) = u[1]*v[2] - u[2]*v[1]

    J = SVector(cross2d(e[1],-e[4]), cross2d(e[2],-e[1]),
                cross2d(e[3],-e[2]), cross2d(e[4],-e[3]))
    maxJ = maximum(abs.(J))
    return maxJ ≈ 0 ? 0.0 : minimum(J) / maxJ
end

################################################################################
### Default Quality Metrics
################################################################################

default_quality_metric(::Simplex) = element_quality_radius_ratio
default_quality_metric(::Block)   = element_quality_mean_ratio

################################################################################
### User-Facing Shorthands
################################################################################

"""
    element_volumes(m::DMesh)

Return a `Vector` of volumes (or areas) for every element in the mesh.
"""
element_volumes(m::DMesh) = collect(element_map(element_volume, m))

"""
    element_qualities(m::DMesh; metric=default_quality_metric(G()))

Return a `Vector` of quality metrics for every element in the mesh. 
"""
function element_qualities(m::DMesh{D, T, G}; metric=default_quality_metric(G())) where {D, T, G}
    return collect(element_map(metric, m))
end

function find_elems(m::DMesh, cond::Function)
    return findall(tt -> cond(sum(m.p[tt]) / length(tt)), m.t)
end


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
function cleanup_mesh(msh::DMesh)
    p, t = msh

    # 1. Snap nodes to a grid to identify duplicates (relative tolerance)
    scaling = maximum(norm.(p))
    if scaling == 0.0
        scaling = 1.0
    end
    pp = [snap.(p1, scaling) for p1 in p]
    
    # 2. Find unique nodes
    ppp = unique(pp)
    
    # 3. Create mappings
    # ix: Mapping from New -> Old (which old node did this new node come from?)
    ix = Int.(indexin(ppp, pp))
    # jx: Mapping from Old -> New (where did this old node go?)
    jx = Int.(indexin(pp, ppp))
    
    # 4. Rebuild Mesh
    new_p = p[ix]
    # Broadcast the index lookup to preserve SVector/Tuple structure of triangles
    new_t = [map(idx -> jx[idx], tri) for tri in t]

    return (msh=DMesh(new_p, new_t), ix=ix)
    
end

"""
    element_face_neighbors(msh::DMesh{D,T,G,N,I}) -> Matrix{Tuple{I, I}}

Compute element connectivities across faces.

Returns a matrix of size `(num_faces_per_element, num_elements)` where each entry 
is a tuple `(neighbor_element, neighbor_local_face)`.

If a face is on the boundary, the entry is `(0, 0)`.
"""
function element_face_neighbors(msh::DMesh{D,T,G,N,I}) where {D,T,G,N,I}
    t = msh.t
    nt = length(t)
    
    fmap = facemap(G())
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
    face_element_map(msh::DMesh{D,T,G,N,I}) -> Dict{SVector{nfv, I}, Vector{Tuple{I, I}}}

Compute the mapping from faces to all connected elements. 

Returns a dictionary where each key is a sorted `SVector` representing the face, 
and the value is a vector of tuples `(element_idx, local_face_idx)`.
"""
function face_element_map(msh::DMesh{D,T,G,N,I}) where {D,T,G,N,I}
    t = msh.t
    nt = length(t)
    
    fmap = facemap(G())
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
function foreach_face(f::Function, msh::DMesh{D,T,G,N,I}) where {D,T,G,N,I}
    nb = element_face_neighbors(msh)
    nt = length(msh.t)
    map = facemap(G())
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
function boundary_faces(msh::DMesh{D,T,G,N,I}) where {D,T,G,N,I}
    map = facemap(G())
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
function all_faces(msh::DMesh{D,T,G,N,I}) where {D,T,G,N,I}
    nfv = length(facemap(G())[1])
    
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
function all_edges(msh::DMesh{D,T,G,N,I}) where {D,T,G,N,I}
    emap = edgemap(G())
    
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
