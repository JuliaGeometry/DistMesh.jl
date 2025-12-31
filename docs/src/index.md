# DistMesh.jl

**DistMesh.jl** is a Julia generator for unstructured triangular and tetrahedral meshes. It uses [Signed Distance Functions](https://en.wikipedia.org/wiki/Signed_distance_function) (SDFs) to define geometries, allowing for clean mesh generation of complex shapes using simple mathematical functions.

## Installation

```julia
using Pkg
Pkg.add("DistMesh")

```

## Quick Start (2D)

The primary entry point for 2D meshing is `distmesh2d`.

First, we generate the mesh.

```@example quickstart
using DistMesh, Plots

msh = distmesh2d(dcircle, huniform, 0.15, ((-1,-1), (1,1)))
```

Now we can visualize the result using `plot`.

```@example quickstart
plot(msh)
```

Finally, if you need the raw coordinate and topology arrays for a solver, you can access the fields directly.

```@example quickstart
# Access the points and connectivity as matrices
p,t = as_arrays(msh)
p
```

```@example quickstart
t
```
