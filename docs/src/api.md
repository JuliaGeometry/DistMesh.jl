# API Reference

## Core Meshing & Types

The primary interface for generating meshes and handling the resulting data structures.

```@docs
distmesh2d
distmesh
DMesh
as_arrays

```

## Distance Functions: Basic Shapes

Pre-defined signed distance functions for common primitives in 2D and 3D.

```@docs
dcircle
drectangle
dhypersphere
dsphere
dblock

```

## Distance Functions: Polygons & Lines

Utilities for working with polygonal boundaries and line segments.

```@docs
dpoly
dline
DistMesh.inpolygon

```

## Distance Functions: Boolean Operations (CSG)

Constructive Solid Geometry operations to combine multiple distance functions.

```@docs
ddiff
dunion
dintersect

```

## Distance Functions: Special Functions

Helper functions for element sizing and specific test cases (like airfoils).

```@docs
dnaca

```

## Mesh utilities: Size functions

```@docs
huniform

```

## Mesh utilities: General

```@docs
element_volumes
element_qualities
cleanup_mesh

```
