module DistMesh

using StaticArrays
using LinearAlgebra
using Delaunator

# --- Load Types ---
include("dmesh.jl")

# --- Load Utilities and 2D Implementation ---
include("distfuncs.jl")
include("meshutils.jl")
include("distmesh2d.jl")

# --- Exports ---

export DMesh, as_arrays
export distmesh2d

export dhypersphere, dcircle, dsphere, drectangle, dblock
export dline, dsegment, dpoly
export ddiff, dunion, dintersect
export huniform
export naca_coeffs, dnaca

export element_qualities, element_volumes, cleanup_mesh

end # module
