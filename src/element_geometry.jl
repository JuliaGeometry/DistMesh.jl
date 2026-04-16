###########################################################################
## Element geometry types

"""Abstract base type for `D`-dimensional element geometries."""
abstract type ElementGeometry{D} end

"""Simplex element (point, line, triangle, tetrahedron)."""
struct Simplex{D} <: ElementGeometry{D} end

"""Tensor-product element (point, line, quadrilateral, hexahedron)."""
struct Block{D} <: ElementGeometry{D} end

const geometry_types = (Simplex, Block)

###########################################################################
## Basic properties

dim(::ElementGeometry{D}) where {D} = D::Int

name(::Simplex{D}) where {D} = D > 3 ? "$(D)D simplex" :
    ("point", "line", "triangle", "tetrahedron")[D+1]
name(::Block{D}) where {D} = D > 3 ? "$(D)D block" :
    ("point", "line", "quadrilateral", "hexahedron")[D+1]

Base.show(io::IO, eg::ElementGeometry) = print(io, "ElementGeometry: $(dim(eg))D $(name(eg))")

nvertices(::ElementGeometry) = error("Not implemented")
nvertices(::Simplex{D}) where {D} = D + 1
nvertices(::Block{D}) where {D} = 2^D

nfaces(::ElementGeometry) = error("Not implemented")
nfaces(::Simplex{D}) where {D} = D + 1
nfaces(::Block{D}) where {D} = 2 * D

nedges(::ElementGeometry) = error("Not implemented")
nedges(::Simplex{D}) where {D} = binomial(D + 1, 2)
nedges(::Block{D}) where {D} = D * 2^(D - 1)

###########################################################################
## Connectivity maps
#
# Maps return an SVector of SVectors.
# Ordering follows standard VTK/ExodusII conventions (CCW/RHR).
# 1D simplices break the opposite-node rule to match 1D blocks.

facemap(::ElementGeometry) = error("Not implemented")
facemap(::Simplex{1}) = SA[SA[1], SA[2]]
facemap(::Simplex{2}) = SA[SA[2, 3], SA[3, 1], SA[1, 2]]
facemap(::Simplex{3}) = SA[SA[2, 3, 4], SA[1, 4, 3], SA[4, 1, 2], SA[3, 2, 1]]

facemap(::Block{1})   = SA[SA[1], SA[2]]
facemap(::Block{2})   = SA[SA[1, 2], SA[2, 3], SA[3, 4], SA[4, 1]]
facemap(::Block{3})   = SA[
    SA[1, 4, 3, 2], # Bottom
    SA[1, 2, 6, 5], # Front
    SA[2, 3, 7, 6], # Right
    SA[3, 4, 8, 7], # Back
    SA[4, 1, 5, 8], # Left
    SA[5, 6, 7, 8]  # Top
]

edgemap(::ElementGeometry) = error("Not implemented")
edgemap(::Simplex{1}) = SA[SA[1, 2]]
edgemap(::Simplex{2}) = SA[SA[2, 3], SA[3, 1], SA[1, 2]]
edgemap(::Simplex{3}) = SA[
    SA[1, 2], SA[2, 3], SA[3, 1], # Base
    SA[1, 4], SA[2, 4], SA[3, 4]  # Pillars to apex
]

edgemap(::Block{1})   = SA[SA[1, 2]]
edgemap(::Block{2})   = SA[SA[1, 2], SA[2, 3], SA[3, 4], SA[4, 1]]
edgemap(::Block{3})   = SA[
    SA[1, 2], SA[2, 3], SA[3, 4], SA[4, 1], # Bottom
    SA[5, 6], SA[6, 7], SA[7, 8], SA[8, 5], # Top
    SA[1, 5], SA[2, 6], SA[3, 7], SA[4, 8]  # Pillars
]

###########################################################################
## Sub-geometry and utilities

"""Return the `newD`-dimensional sub-geometry type of `eg`."""
subgeom(::Simplex, newD) = Simplex{newD}()
subgeom(::Block,   newD) = Block{newD}()

"""Infer element geometry from spatial dimension `D` and vertex count `nv`."""
function find_elgeom(D, nv)
    nv == nvertices(Block{D}())   && return Block{D}()
    nv == nvertices(Simplex{D}()) && return Simplex{D}()
    nv == nvertices(Block{D-1}())   && return Block{D-1}()
    nv == nvertices(Simplex{D-1}()) && return Simplex{D-1}()
    error("Cannot determine element geometry for D=$D, nv=$nv")
end
