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

function simplex_area(el)
    if length(el) == 3 # Triangle
        p12 = el[2] - el[1]
        p13 = el[3] - el[1]
        return (p12[1] * p13[2] - p12[2] * p13[1]) / 2
    else
        @assert "Dimension not implemented"
    end
end

function simplex_qual(el)
    if length(el) == 3 # Triangle
        norm(vec) = sqrt(sum(vec.^2))
        a,b,c = ( norm(el[ix[2]] - el[ix[1]]) for ix in ((1,2),(2,3),(3,1)) )
        r = 0.5*sqrt((b+c-a) * (c+a-b) * (a+b-c) / (a+b+c))
        R = a*b*c / sqrt((a+b+c) * (b+c-a) * (c+a-b) * (a+b-c))
        return 2*r/R
    else
        @assert "Dimension not implemented"
    end
end

elemwise_feval(m::DMesh, f) = [ f(m.p[tt]) for tt in m.t ]

simpqual(m::DMesh, fqual=simplex_qual) = elemwise_feval(m, fqual)
simpvol(m::DMesh) = elemwise_feval(m, simplex_area)

################################################################################
### General mesh utilities
################################################################################

snap(x::T, scaling=1) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling*tol*round(x/scaling/tol) + zero(T)  # Adding zero to uniquify -0.0 and 0.0


"""
    fixmesh(msh::DMesh) -> (msh::DMesh, ix::Vector{Int})

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
clean_msh, = fixmesh(dirty_msh)  # Ignoring the index output (ix)

```

"""
function fixmesh(msh::DMesh)
    p, t = msh

    # 1. Snap nodes to a grid to identify duplicates (relative tolerance)
    scaling = maximum(norm.(p))
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
