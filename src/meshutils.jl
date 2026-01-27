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
