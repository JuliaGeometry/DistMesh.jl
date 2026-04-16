module DistMesh

using StaticArrays
using LinearAlgebra
using Delaunator

# --- Load Types ---
include("element_geometry.jl")
include("dmesh.jl")

# --- Load Utilities and 2D Implementation ---
include("distfuncs.jl")
include("meshutils.jl")
include("io.jl")
include("plotting.jl")
include("distmesh2d.jl")

# --- Exports ---

export ElementGeometry, Simplex, Block
export nvertices, nfaces, nedges, facemap, edgemap

export DMesh, as_arrays
export distmesh2d
export get_camera_view

export dhypersphere, dcircle, dsphere, drectangle, dblock
export dline, dsegment, dpoly
export ddiff, dunion, dintersect
export huniform
export naca_coeffs, dnaca

export read_stl

export element_qualities, element_volumes, find_elems, cleanup_mesh
export element_face_neighbors, face_element_map, find_boundary_elements, find_nonmanifold_elements
export all_faces, boundary_faces, boundary_nodes, all_edges, node_degrees

end # module
