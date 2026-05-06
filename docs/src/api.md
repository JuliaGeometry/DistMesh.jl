# API Reference

## Core Meshing & Types

```@docs
distmesh2d
DMesh
as_arrays
```

## Distance Functions: Basic Shapes

```@docs
dcircle
drectangle
dhypersphere
dsphere
dblock
```

## Distance Functions: Polygons & Lines

```@docs
dpoly
dline
DistMesh.inpolygon
```

## Distance Functions: Boolean Operations (CSG)

```@docs
ddiff
dunion
dintersect
```

## Distance Functions: Special Functions

```@docs
dnaca
```

## Mesh Utilities: Size Functions

```@docs
huniform
```

## Mesh Utilities: General

```@docs
element_volumes
element_qualities
cleanup_mesh
is_manifold_mesh
```

## Element Utilities

```@docs
element_face_neighbors
face_element_map
find_boundary_elements
find_nonmanifold_elements
```

## Face and Edge Utilities

```@docs
all_faces
boundary_faces
all_edges
```

## Node Utilities

```@docs
boundary_nodes
node_degrees
node_adjacency
node_element_map
```
