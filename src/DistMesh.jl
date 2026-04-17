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
include("quality_metrics.jl")
include("topology.jl")
include("mesh_improvement.jl")
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

export read_stl, write_stl, read_ply, write_ply

export element_qualities, element_volumes, find_elems, cleanup_mesh
export element_face_neighbors, face_element_map, find_boundary_elements, find_nonmanifold_elements, is_manifold_mesh
export all_faces, boundary_faces, boundary_nodes, all_edges, node_degrees
export trimesh_flip!, trimesh_collapse

end # module
