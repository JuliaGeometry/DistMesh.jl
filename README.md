# DistMesh.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://distmesh.juliageometry.org/dev)
[![Codecov](https://codecov.io/gh/juliageometry/DistMesh.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliageometry/DistMesh.jl)

**DistMesh.jl** is a Julia package for generating unstructured triangular and tetrahedral meshes using Signed Distance Functions (SDFs). 

It is designed for simplicity and interactive visualization, making it easy to generate high-quality meshes for finite element analysis or geometric modeling.

## Quick Start (2D)

The core function `distmesh2d` generates a mesh based on a distance function (negative inside, positive outside).

```julia
using DistMesh
using GLMakie # Optional, for plotting

# 1. Define the geometry
# (Here we use the built-in 'dcircle' and 'huniform' size function)
# Region: Box from (-1,-1) to (1,1)
# Initial edge length: 0.2
msh = distmesh2d(dcircle, huniform, 0.2, ((-1,-1), (1,1)))

# 2. Visualize
# DistMesh automatically plots correctly with Makie or Plots.jl
plot(msh)

```

### Key Features

* **Simple API:** Define geometry using standard Julia functions.
* **Robust:** Uses a physical force-equilibrium model to optimize node locations.
* **Visual:** Built-in `show` methods and plotting recipes for instant feedback.

---

## Legacy & 3D Support

### N-Dimensional Generation

The package includes the original N-dimensional implementation (supporting 3D and higher) via the internal `DistMeshND` module. This implementation closely mirrors the original MATLAB source.

### Background

This package is a Julia port of the [DistMesh](http://persson.berkeley.edu/distmesh/) algorithm developed by [Per-Olof Persson](http://persson.berkeley.edu/). Several improvements have been made to performance and type stability compared to the original version.

[Technical Report](https://sjkellyorg.files.wordpress.com/2020/11/distmesh_sjkelly.pdf)
