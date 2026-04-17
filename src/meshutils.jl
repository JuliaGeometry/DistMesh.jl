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
### User-Facing Shorthands
################################################################################

"""
    element_volumes(m::DMesh)

Return a `Vector` of volumes (or areas) for every element in the mesh.
"""
element_volumes(m::DMesh) = collect(element_map(element_volume, m))

function find_elems(m::DMesh, cond::Function)
    return findall(tt -> cond(sum(m.p[tt]) / length(tt)), m.t)
end

