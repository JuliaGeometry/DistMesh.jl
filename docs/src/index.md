# DistMesh.jl

**DistMesh.jl** is a Julia generator for unstructured triangular and tetrahedral meshes. It uses [Signed Distance Functions](https://en.wikipedia.org/wiki/Signed_distance_function) (SDFs) to define geometries, enabling the generation of high-quality, isotropic meshes for complex shapes defined by simple mathematical functions.

## Installation

```julia
using Pkg
Pkg.add("DistMesh")
```

## Introduction

### A Simple Example Mesh

The primary entry point for 2D meshing is `distmesh2d`. It generates a mesh based on:

* A distance function (negative inside, positive outside). In this example,
  we represent the unit circle by $d(x,y) = \sqrt{x^2+y^2} - 1$
  
* A relative size function. For a uniform mesh, we can simply set $h(x,y) = 1$

* An initial, or minimum, element edge length. Since our size function is
  uniform, this will also be roughly the size of the generated elements.
  
* A bounding box for the domain. For the unit circle, we use the coordinates
  (-1,-1) and (1,1).
  
The code below demonstrates how to implement this using DistMesh in Julia.

```@example introduction
using DistMesh
using GLMakie                # or Plots, or CairoMakie (optional)

fd(p) = sqrt(sum(p.^2)) - 1  # or dcircle(p) - unit circle geometry
fh(p) = 1.0                  # or huniform(p) - uniform size function
hmin  = 0.2                  # initial edge lengths
bbox  = ((-1,-1), (1,1))     # bounding box for unit circle

msh = distmesh2d(fd, fh, hmin, bbox)
```

The resulting mesh can be visualized with the `plot` command.

```@example introduction
plot(msh)
```

### Accessing the generated mesh

The `msh` object contains the node coordinates `p` and the triangle indices `t`, stored as vectors of static vectors (for efficiency)

```@example introduction
# Access the nodes and connectivity as vectors of static vectors
p,t = msh
p[1:3] # First 3 nodes
```

```@example introduction
t[1:3] # First 3 triangles
```

If you prefer a matrix-based representation (similar to the original MATLAB DistMesh, except transposed), the utility function `as_arrays` creates a zero-allocation view for this. Note that if you modify these matrices, the original versions in `msh` will also be changed!

```@example introduction
# Access the nodes and connectivity as matrices
p_mat,t_mat = as_arrays(msh)
p_mat[:,1:3] # First 3 nodes
```

```@example introduction
t_mat[:,1:3] # First 3 triangles
```

If you want separate re-allocated versions of these matrices, use `collect`.

### Fixed points

An optional argument to `distmesh2d` is `pfix` - a vector of frozen points that are part of the generated mesh. These are typically required for boundaries with sharp corners, since the implicit geometry representation as a signed distance function does not explicitly provide them. Here is a simple example of a mesh of a rectangle, where we include the four corner points in `pfix1`.
```@example introduction
fd(p) = drectangle(p, -1, 1, -1, 1)
pfix = ((-1,-1), (1,-1), (-1,1), (1,1))
# fh, hmin, bbox same as before
msh = distmesh2d(fd, fh, hmin, bbox, pfix)
```

```@example introduction
plot(msh)
```

## Distance Functions

## Size Functions

