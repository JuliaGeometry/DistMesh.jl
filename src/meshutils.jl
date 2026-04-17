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

using LinearAlgebra

################################################################################
### Utilities
################################################################################

# Helper to compute the magnitude of the cross product for both 2D and 3D SVectors.
_cross_mag(u::SVector{2}, v::SVector{2}) = abs(u[1]*v[2] - u[2]*v[1])
_cross_mag(u::SVector{3}, v::SVector{3}) = norm(cross(u, v))

################################################################################
### Element Volumes
################################################################################

element_volume(::ElementGeometry, el) = error("Not implemented for this geometry")

element_volume(::Simplex{1}, el) = norm(el[2] - el[1])
element_volume(::Block{1}, el)   = norm(el[2] - el[1])

"""
    element_volume(::Simplex{2}, el)

Compute the area of a 2D triangle or 3D surface triangle.
"""
function element_volume(::Simplex{2}, el)
    p1, p2, p3 = el
    return _cross_mag(p2 - p1, p3 - p1) / 2
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

Compute the area of a 2D quadrilateral or 3D surface quadrilateral.
(Calculated as half the magnitude of the cross product of its diagonals).
"""
function element_volume(::Block{2}, el)
    p1, p2, p3, p4 = el
    return _cross_mag(p3 - p1, p4 - p2) / 2
end

"""
    element_volume(::Block{3}, el)

Compute the volume of a 3D hexahedron.
"""
function element_volume(::Block{3}, el)
    p1, p2, p3, p4, p5, p6, p7, p8 = el
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
    a, b, c = norm(p2 - p1), norm(p3 - p2), norm(p1 - p3)
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

function element_quality_mean_ratio(::Block{2}, el)
    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4)

    l2 = SVector(sum(abs2, e[1]), sum(abs2, e[2]), sum(abs2, e[3]), sum(abs2, e[4]))
    A  = SVector(_cross_mag(e[1], -e[4]), _cross_mag(e[2], -e[1]),
                 _cross_mag(e[3], -e[2]), _cross_mag(e[4], -e[3]))
                 
    Q  = SVector((l2[1]+l2[3]), (l2[2]+l2[4]),
                 (l2[3]+l2[1]), (l2[4]+l2[2])) ./ (2 .* A)

    return 4 / sum(Q)
end

function element_quality_condition_number(::Block{2}, el)
    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4) 

    l2   = SVector(sum(abs2, e[1]), sum(abs2, e[2]), sum(abs2, e[3]), sum(abs2, e[4]))
    sins = SVector(_cross_mag(e[1], -e[4]), _cross_mag(e[2], -e[1]),
                   _cross_mag(e[3], -e[2]), _cross_mag(e[4], -e[3])) ./
           SVector(sqrt(l2[4]*l2[1]), sqrt(l2[1]*l2[2]),
                   sqrt(l2[2]*l2[3]), sqrt(l2[3]*l2[4]))
                   
    any(<=(0), sins) && return 0.0

    k = SVector(l2[4]+l2[1], l2[1]+l2[2], l2[2]+l2[3], l2[3]+l2[4]) ./
        (SVector(sqrt(l2[4]*l2[1]), sqrt(l2[1]*l2[2]),
                 sqrt(l2[2]*l2[3]), sqrt(l2[3]*l2[4])) .* sins)

    return 4 / sqrt(sum(abs2, k))
end

function element_quality_min_scaled_jacobian(::Block{2}, el)
    length(first(el)) == 2 || error("Minimum scaled jacobian for Block{2} is strictly defined for 2D meshes. Surface quad implementation (3D) is not supported.")

    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4)  
    cross2d(u, v) = u[1]*v[2] - u[2]*v[1]

    J = SVector(cross2d(e[1],-e[4]), cross2d(e[2],-e[1]),
                cross2d(e[3],-e[2]), cross2d(e[4],-e[3]))
    maxJ = maximum(abs.(J))
    return maxJ ≈ 0 ? 0.0 : minimum(J) / maxJ
end

function element_quality_mean_ratio(::Block{3}, el)
    p1,p2,p3,p4,p5,p6,p7,p8 = el
    corners = (
        (p2-p1, p4-p1, p5-p1), (p3-p2, p1-p2, p6-p2),
        (p4-p3, p2-p3, p7-p3), (p1-p4, p3-p4, p8-p4),
        (p8-p5, p6-p5, p1-p5), (p5-p6, p7-p6, p2-p6),
        (p6-p7, p8-p7, p3-p7), (p7-p8, p5-p8, p4-p8),
    )
    function corner_quality(e1, e2, e3)
        detW = dot(e1, cross(e2, e3))
        detW <= 0 && return 0.0
        frob2 = sum(abs2, e1) + sum(abs2, e2) + sum(abs2, e3)
        return 3 * cbrt(detW^2) / frob2
    end
    return minimum(corner_quality(c...) for c in corners)
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
function cleanup_mesh(msh::DMesh{D,T,G,N,I}) where {D,T,G,N,I}
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


"""
    trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}; tol=1e-3) -> DMesh
    trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}, nb; tol=1e-3) -> DMesh

Collapse short edges in a triangular mesh, returning a new mesh with fewer vertices and triangles.

Any edge with Euclidean length below `tol` is collapsed: one endpoint is removed and all
references to it are redirected to the surviving endpoint. When one endpoint is a boundary
vertex and the other is not, the boundary vertex is always kept. Otherwise the lower-index
vertex is kept.

For interior edges, both adjacent triangles are removed. For boundary edges, only the one
triangle containing the edge is removed. Neighboring triangles that become degenerate
(two or more vertices mapped to the same point) are also filtered out.

A pre-computed neighbor matrix `nb` from `element_face_neighbors` may be supplied to
avoid recomputing it.

# Arguments
- `msh`: Input triangular mesh (`Simplex{2}` elements, any spatial dimension)
- `nb`: Optional pre-computed element face neighbor matrix
- `tol`: Edge length threshold below which edges are collapsed (default: `1e-3`)

# Returns
- A new `DMesh` with short edges removed.

# Example
```julia
msh2 = trimesh_collapse(msh, tol=0.01)
```
"""
function trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}; tol=1e-3) where {D,T}
    nb = element_face_neighbors(msh)
    trimesh_collapse(msh, nb; tol)
end

function trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}, nb; tol=1e-3) where {D,T}
    np = length(msh.p)
    nt = length(msh.t)
    fmap = facemap(Simplex{2}())

    # --- Phase 1: Identify boundary vertices ---
    is_boundary = falses(np)
    for iel in 1:nt
        for jf in 1:3
            if nb[jf, iel][1] == 0
                edge = msh.t[iel][fmap[jf]]
                is_boundary[edge[1]] = true
                is_boundary[edge[2]] = true
            end
        end
    end

    # --- Phase 2: Mark collapses ---
    remap = collect(1:np)           # remap[i] = surviving vertex for i
    dead_tris = falses(nt)

    # Helper: follow remap chain to find current live representative
    function resolve(v)
        while remap[v] != v
            v = remap[v]
        end
        return v
    end

    for iel in 1:nt
        for jf in 1:3
            t2, _ = nb[jf, iel]

            # Skip interior edges where we'd process the same edge from t2 side
            t2 > 0 && t2 < iel && continue

            edge = msh.t[iel][fmap[jf]]
            va = resolve(edge[1])
            vb = resolve(edge[2])

            # Already collapsed to same point
            va == vb && continue

            # Check edge length
            norm(msh.p[va] - msh.p[vb]) >= tol && continue

            # Decide which vertex to keep: prefer boundary vertex
            v_keep, v_del = if is_boundary[va] && !is_boundary[vb]
                va, vb
            elseif is_boundary[vb] && !is_boundary[va]
                vb, va
            else
                min(va, vb), max(va, vb)
            end

            # Apply collapse
            remap[v_del] = v_keep

            # Mark triangles as dead
            dead_tris[iel] = true
            if t2 > 0
                dead_tris[t2] = true
            end
        end
    end

    # --- Phase 3: Path-compress remap chains (A→B→C becomes A→C) ---
    for i in 1:np
        remap[i] = resolve(i)
    end

    # --- Phase 4: Build compacted vertex list ---
    # Surviving vertices: those that map to themselves
    keep_verts = findall(i -> remap[i] == i, 1:np)

    # new_idx[old_index] = new index in compacted array (only valid for kept verts)
    new_idx = zeros(Int, np)
    for (k, v) in enumerate(keep_verts)
        new_idx[v] = k
    end

    # Composed map: old vertex i → new index
    final_idx(i) = new_idx[remap[i]]

    # --- Phase 5: Build compacted element list ---
    new_t = SVector{3, Int}[]
    sizehint!(new_t, nt - count(dead_tris))

    for iel in 1:nt
        dead_tris[iel] && continue

        tri = msh.t[iel]
        a, b, c = final_idx(tri[1]), final_idx(tri[2]), final_idx(tri[3])

        # Drop degenerate elements (two or more vertices collapsed to same point)
        (a == b || b == c || a == c) && continue

        push!(new_t, SVector(a, b, c))
    end

    return DMesh(msh.p[keep_verts], new_t)
end

"""
    trimesh_flip!(msh::DMesh{D,T,Simplex{2},3})
    trimesh_flip!(msh::DMesh{D,T,Simplex{2},3}, nb)

Improve triangular mesh quality by flipping edges between adjacent triangles.

Modifies the mesh in-place by examining each edge shared between two adjacent triangles
and flipping it across the opposite diagonal if the quality of both triangles improves
significantly. This operation is commonly used in mesh optimization algorithms.

The quality of a flip is evaluated using the mean ratio quality metric. An edge is flipped
if the minimum quality of the two new triangles exceeds the minimum quality of the original
pair by at least 0.025. For 3D meshes, the flip is additionally checked to ensure that
the resulting triangles have consistent normal orientations.

A pre-computed neighbor matrix `nb` from `element_face_neighbors` may be supplied to
avoid recomputing it.

# Arguments
- `msh`: Input triangular mesh (`Simplex{2}` elements, 2D or 3D)
- `nb`: Optional pre-computed element face neighbor matrix

# Example
```julia
trimesh_flip!(msh)  # Improve mesh quality in-place
```
"""
function trimesh_flip!(msh::DMesh{D,T,Simplex{2},3}) where {D,T}
    nb = element_face_neighbors(msh)
    trimesh_flip!(msh, nb)
end

function trimesh_flip!(msh::DMesh{D,T,Simplex{2},3}, nb) where {D,T}
    function trinormal3(pts)
        p1, p2, p3 = pts
        c = cross(p2 - p1, p3 - p1)
        return c / norm(c)
    end
    
    triqual3(pts) = element_quality_mean_ratio(Simplex{2}(), pts)

    nt = length(msh.t)
    for t1 in 1:nt
        for n1 in 1:3
            t2, n2 = nb[n1, t1]
            
            if t2 > 0
                old_t1 = msh.t[t1]
                old_t2 = msh.t[t2]
                
                q1 = triqual3(msh.p[old_t1])
                q2 = triqual3(msh.p[old_t2])
                minqold = min(q1, q2)
                
                if minqold < 0.9
                    tix11 = mod(n1, 3) + 1
                    tix12 = mod(n1 + 1, 3) + 1
                    tix21 = mod(n2, 3) + 1
                    tix22 = mod(n2 + 1, 3) + 1

                    # Create new triangles with swapped edges (no allocations with ntuple)
                    newt1 = SVector(ntuple(i -> i == tix12 ? old_t2[n2] : old_t1[i], 3))
                    newt2 = SVector(ntuple(i -> i == tix22 ? old_t1[n1] : old_t2[i], 3))

                    q3 = triqual3(msh.p[newt1])
                    q4 = triqual3(msh.p[newt2])
                    minqnew = min(q3, q4)

                    if minqnew > minqold + 0.025
                        flip = if D == 2
                            true  # In 2D, quality function handles orientation
                        else  # D == 3
                            normal1 = trinormal3(msh.p[old_t1])
                            normal2 = trinormal3(msh.p[old_t2])
                            normal3 = trinormal3(msh.p[newt1])
                            normal4 = trinormal3(msh.p[newt2])
                            if minqold < 0.001 # Desperate
                                (dot(normal3, normal4) > 0)
                            else
                                (dot(normal1, normal2) > 0) && (dot(normal3, normal4) > 0)
                            end
                        end
                        
                        if flip
                            # Update triangles in mesh
                            msh.t[t1] = newt1
                            msh.t[t2] = newt2

                            # Update neighbor connectivity
                            nbt, nbn = nb[tix21, t2]
                            nb[n1, t1] = (nbt, nbn)
                            if nbt > 0
                                nb[nbn, nbt] = (t1, n1)
                            end

                            nbt, nbn = nb[tix11, t1]
                            nb[n2, t2] = (nbt, nbn)
                            if nbt > 0
                                nb[nbn, nbt] = (t2, n2)
                            end

                            nb[tix11, t1] = (t2, tix21)
                            nb[tix21, t2] = (t1, tix11)
                        end
                    end
                end
            end
        end
    end
end
