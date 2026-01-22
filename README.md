# DistMesh.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://distmesh.juliageometry.org/dev)
[![Codecov](https://codecov.io/gh/juliageometry/DistMesh.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliageometry/DistMesh.jl)

**DistMesh.jl** is a Julia generator for unstructured triangular and tetrahedral meshes. It uses [Signed Distance Functions](https://en.wikipedia.org/wiki/Signed_distance_function) (SDFs) to define geometries, enabling the generation of high-quality, isotropic meshes for complex shapes defined by simple mathematical functions.

Primary use cases include Finite Element Analysis (FEA), computational fluid dynamics, and geometric modeling.

---

## Note on Version 0.2.0

Version 0.2.0 represents a major rewrite of `DistMesh.jl`, focusing on replicating the original 2D code but with better performance, native Julia implementations, and integration with the Makie/Plots ecosystems. **This version is now fully MIT licensed.**

**Legacy Support:** The original N-dimensional meshing code port has been retired in this version to allow for a cleaner architecture. If you require the legacy `distmeshnd` functionality, please pin your package version to v0.1:

```julia
] add DistMesh@0.1
```

---

## Quick Start

The core function `distmesh2d` generates a mesh based on a distance function, a relative size function, an initial edge length, and a bounding box.

The code below generates a mesh of the unit circle:

```julia
using DistMesh

# Generate mesh: Distance function, size function, resolution, bounding box
msh = distmesh2d(dcircle, huniform, 0.2, ((-1,-1), (1,1)))
```

```text
DMesh: 2D, 88 nodes, 143 triangle elements
```

Optionally, the mesh can be visualized using various plotting packages:

```julia
using GLMakie # or Plots, or CairoMakie

plot(msh)
```

![Triangular mesh of the unit circle](circle_mesh.png)

For more details and extensive examples, please see the [DistMesh documentation](https://distmesh.juliageometry.org).

---

## Background

This package is a Julia port of the [DistMesh](http://persson.berkeley.edu/distmesh/) algorithm developed by [Per-Olof Persson](http://persson.berkeley.edu/). Significant improvements have been made to performance and type stability compared to the original MATLAB implementation. The algorithm is described in the following publications:

* P.-O. Persson, G. Strang, *[A Simple Mesh Generator in MATLAB](https://persson.berkeley.edu/distmesh/persson04mesh.pdf)*. SIAM Review, Volume 46 (2), pp. 329-345, June 2004.
* P.-O. Persson, *[Mesh Generation for Implicit Geometries](https://persson.berkeley.edu/thesis/persson-thesis-color.pdf)*. Ph.D. thesis, Department of Mathematics, MIT, Dec 2004.

---

## Related Packages

Several other implementations of the DistMesh algorithm exist in the Julia ecosystem, including [DistMesh-Julia](https://github.com/precise-simulation/distmesh-julia) and [DistMesh2D.jl](https://juliapackages.com/p/distmesh2d). Please consider checking these packages if DistMesh.jl does not meet your specific requirements.
