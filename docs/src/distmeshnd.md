# Legacy & Theory

The `DistMeshND` module preserves the original N-dimensional algorithm structure, closely mirroring the original MATLAB implementation.

## Background & Theory

DistMesh.jl implements **simplex refinement** on signed distance functions. The algorithm was first presented in 2004 by Per-Olof Persson and assumes the geometry is defined by any function that returns a signed distance (negative inside, positive outside).

### What is Simplex Refinement?

In layman's terms, a **simplex** is a triangle in 2D or a tetrahedron in 3D. When simulating physics, you often want a mesh of simplices that offers:

* **Accurate approximation** of boundaries and features.
* **Adaptive mesh sizes** to improve accuracy where needed.
* **High Element Quality** (Near-regular simplices).

DistMesh is designed to address these needs using a physical analogy: it treats the mesh nodes as particles connected by springs (edges). The nodes "relax" into a configuration that balances the internal spring forces with the geometric constraints of the Signed Distance Function.

### Algorithm Overview

1. **Initial Distribution:** Nodes are distributed stochastically or on a grid inside the bounding box.
2. **Delaunay Triangulation:** Nodes are connected to form a mesh topology.
3. **Force Equilibrium:** Edges act as springs. Nodes move to minimize energy.
4. **Boundary Projection:** Nodes that move outside the geometry are projected back onto the zero-level set (the boundary) using the SDF gradient.
5. **Refinement:** Triangles that are too large (according to your size function) are split; edges that are too short are collapsed.

### Comparison to Other Methods

DistMesh generally has a very low memory footprint and avoids complex boundary patching logic found in advancing front methods. Since the global state of simplex qualities is optimized in each iteration, this often leads to very high-quality meshes for smooth geometries.

Unlike surface-based meshing tools (which require Piecewise Linear Complexes/STLs as input), DistMesh works directly on mathematical functions or volumetric data, making it ideal for topology optimization, level-set methods, and image-based meshing.

### Differences from MATLAB

Given the same parameters, the Julia implementation of DistMesh will generally perform significantly faster than the original MATLAB implementation.

* **Memory:** Reduced allocation overhead.
* **Triangulation:** While MATLAB uses QHull, Julia's ecosystem (and this package's legacy N-D modules) often leverages `VoronoiDelaunay.jl` or `MiniQhull` for efficiency.

### Working with Signed Distance Functions

You can define SDFs manually (as shown in the examples), or generate them from data. Useful libraries for turning gridded/level-set data into SDFs include:

* [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
* AdaptiveDistanceFields.jl

---

## Legacy API

These functions belong to the internal `DistMeshND` module.

```@autodocs
Modules = [DistMesh.DistMeshND]
```
