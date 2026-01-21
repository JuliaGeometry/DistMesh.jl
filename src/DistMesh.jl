module DistMesh

using StaticArrays
using LinearAlgebra
using Delaunator

# --- 1. Load Types ---
include("dmesh.jl")

# --- 2. Load Utilities and 2D Implementation ---
include("distfuncs.jl")
include("meshutils.jl")
include("distmesh2d.jl")

# --- 3. Load Legacy 3D Module ---
include("distmeshnd/DistMeshND.jl")
using .DistMeshND

# --- 4. Exports ---

# New API
export DMesh, as_arrays
export distmesh2d

export dhypersphere, dcircle, dsphere, drectangle, dblock
export dline, dsegment, dpoly
export ddiff, dunion, dintersect
export huniform
export naca_coeffs, dnaca

export element_qualities, element_volumes, cleanup_mesh

# Legacy API (Re-exports)
export distmeshnd, DistMeshSetup, DistMeshStatistics, HUniform

# --- 5. Deprecations ---

"""
    distmesh(args...)

[DEPRECATED] Please use `distmeshnd` for 3D meshes or `distmesh2d` for 2D meshes.
"""
function distmesh(args...; kwargs...)
    Base.depwarn("distmesh() is deprecated. Please use `distmeshnd()` or `distmesh2d()`.", :distmesh)
    return distmeshnd(args...; kwargs...)
end

export distmesh

end # module
