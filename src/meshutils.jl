################################################################################
### Basic mesh element properties
################################################################################

const facemap = [SA[SA[2,1]],
                 SA[SA[2,3],
                    SA[3,1],
                    SA[1,2]],
                 SA[SA[2,3,4],
                    SA[1,4,3],
                    SA[4,1,2],
                    SA[3,2,1]]]

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
### Element properties
################################################################################

"""
    element_volume(el)

Compute the generalized volume (area in 2D, volume in 3D) of a single element.
`el` is a vector of coordinates (e.g. [[x1, y1], [x2, y2], [x3, y3]]).
"""
function element_volume(el)
    # Generic simplex handling could go here later.
    # For now, explicit checks for standard shapes:

    N = length(el)
    D = N > 0 ? length(el[1]) : 0
    
    if D == 2 && N == 3 # Triangle (2D)
        p12 = el[2] - el[1]
        p13 = el[3] - el[1]
        return (p12[1] * p13[2] - p12[2] * p13[1]) / 2
        
    elseif D == 3 && N == 4 # Tetrahedron (3D)
        error("3D Tetrahedra not yet implemented")
        
    elseif D == 2 && N == 4 # Example: Quad (2D) logic check
        # Quad logic
        error("2D Quadrilaterals not yet implemented")
        
    else
        error("Element type or dimension not supported")
    end
end

"""
    element_quality(el)

Compute a quality metric for a single element (normalized 0.0 to 1.0).
Currently implements 2*r/R (radius ratio) for triangles.
"""
function element_quality(el)
    if length(el) == 3 # Triangle
        # Lengths of the three edges
        a = norm(el[2] - el[1])
        b = norm(el[3] - el[2])
        c = norm(el[1] - el[3])
        
        # Semiperimeter
        s = (a + b + c) / 2
        
        # Area (Heron's formula) for inradius calculation
        # area = sqrt(s * (s-a) * (s-b) * (s-c))
        # r = area / s
        # R = a*b*c / (4*area)
        # Quality = 2*r/R
        
        denom = (a * b * c)
        if denom ≈ 0
             return 0.0
        end
        
        return 8 * (s - a) * (s - b) * (s - c) / denom
    else
         error("Dimension not implemented")
    end
end

# Helper to map a function over all elements
_map_elements(m::DMesh, f) = [f(m.p[indices]) for indices in m.t]

"""
    element_qualities(m::DMesh, quality_func=element_quality)

Return a vector of quality metrics for every element in the mesh.
"""
element_qualities(m::DMesh, f=element_quality) = _map_elements(m, f)

"""
    element_volumes(m::DMesh)

Return a vector of volumes (or areas) for every element in the mesh.
"""
element_volumes(m::DMesh) = _map_elements(m, element_volume)

################################################################################
### General mesh utilities
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
    element_face_neighbors(msh::DMesh{D,T,N,I}) -> (t2t, t2n)

Compute element connectivities (topology) across faces.

Returns two matrices of size `(num_faces_per_element, num_elements)`:
- `t2t[j, i]`: The index of the neighbor element sharing face `j` of element `i`.
- `t2n[j, i]`: The local face index of that neighbor element (seen from the neighbor).

If element `i`'s face `j` is on the boundary, both entries are `0`.

# Arguments
- `msh`: The mesh object.

# Returns
- `t2t`: Element-to-neighbor-element map.
- `t2n`: Element-to-neighbor-face map.
"""
function element_face_neighbors(msh::DMesh{D,T,N,I}) where {D,T,N,I}
    t = msh.t
    nt = length(t)
    
    map = facemap[N-1] 
    nf = length(map)       # Number of faces per element
    nfv = length(map[1])   # Number of vertices per face

    t2t = zeros(I, nf, nt)
    t2n = zeros(I, nf, nt)
    
    # Key is SVector (canonical face), Value is (element_idx, face_idx)
    dd = Dict{SVector{nfv, I}, NTuple{2, I}}()
    sizehint!(dd, nt * nf)

    for iel in 1:nt
        verts = t[iel]
        
        for jf in 1:nf
            # Construct the canonical face key
            key = sort(verts[map[jf]])
            
            if haskey(dd, key)
                nbel, nbface = pop!(dd, key)
                
                # Link Current -> Neighbor
                t2t[jf, iel] = nbel
                t2n[jf, iel] = nbface
                
                # Link Neighbor -> Current
                t2t[nbface, nbel] = iel
                t2n[nbface, nbel] = jf
            else
                dd[key] = (iel, jf)
            end
        end
    end
    return t2t, t2n
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
function boundary_faces(msh::DMesh{D,T,N,I}) where {D,T,N,I}
    t2t, = element_face_neighbors(msh)
    nt = length(msh.t)
    
    map = facemap[N-1] 
    nf = length(map)       # Number of faces per element
    nfv = length(map[1])   # Number of vertices per face (L in docstring)

    bnd = SVector{nfv,I}[]
    sizehint!(bnd, floor(Int, sqrt(nt) * 4)) 

    for iel in 1:nt
        for jf in 1:nf
            if t2t[jf,iel] == 0
                push!(bnd, msh.t[iel][map[jf]])
            end
        end
    end
    return bnd
end
