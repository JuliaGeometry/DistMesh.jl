# DistMesh.jl

**DistMesh.jl** is a Julia package for generating unstructured triangular and tetrahedral meshes. It uses [Signed Distance Functions](https://en.wikipedia.org/wiki/Signed_distance_function) (SDFs) to define geometries, enabling the generation of high-quality, isotropic meshes for complex shapes defined by simple mathematical functions.

---

## Installation

```julia
using Pkg
Pkg.add("DistMesh")

```

---

## Introduction

### A Simple Example Mesh

The primary entry point for 2D meshing is `distmesh2d`. It generates a mesh based on:

* A **distance function** `fd(p)` (negative inside, positive outside). In this example, we represent the unit circle by $d(x,y) = \sqrt{x^2+y^2} - 1$.
* A **relative size function** `fh(p)`. For a uniform mesh, we can simply set $h(x,y)=1$.
* An **initial element edge length** `h0`. Since our size function is uniform, this will also be roughly the final size of the generated elements.
* A **bounding box** `bbox` for the domain. For the unit circle, we use the coordinates `((-1,-1), (1,1))`.

The code below demonstrates how to implement this using DistMesh in Julia.

```@example introduction
using DistMesh
using CairoMakie             # or Plots, or GLMakie (optional)
CairoMakie.activate!(type="png", px_per_unit=1.0) # hide

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

The `msh` object contains the node coordinates `p` and the triangle indices `t`, stored as vectors of static vectors (for efficiency).

```@example introduction
# Access the nodes and connectivity as vectors of static vectors
p, t = msh
p[1:3] # First 3 nodes

```

```@example introduction
t[1:3] # First 3 triangles (indices)

```

```@example introduction
p[t[1]] # x,y-coordinates of the first triangle

```

If you prefer a matrix-based representation (similar to the original MATLAB DistMesh, except transposed), the utility function `as_arrays` creates a zero-allocation view for this. Note that if you modify these matrices, the original versions in `msh` will also be changed!

```@example introduction
# Access the nodes and connectivity as matrices
p_mat, t_mat = as_arrays(msh)
p_mat[:, 1:3] # First 3 nodes

```

```@example introduction
t_mat[:, 1:3] # First 3 triangles (indices)

```

```@example introduction
p_mat[:, t_mat[:, 1]] # x,y-coordinates of the first triangle

```

If you want separate re-allocated versions of these matrices, use `collect`.

### Fixed points

An optional argument to `distmesh2d` is `pfix` - a vector of frozen points that are forced to be part of the generated mesh. These are typically required for boundaries with sharp corners, since the implicit geometry representation as a signed distance function does not explicitly provide them.

Here is a simple example of a mesh of a rectangle, where we include the four corner points in `pfix`.

```@example introduction
fd(p) = drectangle(p, 0, 1, 0, 1)
pfix = ((0,0), (1,0), (0,1), (1,1))
msh = distmesh2d(fd, huniform, 0.1, ((0,0), (1,1)), pfix)
plot(msh)

```

### Implicit Geometry vs Distance Function

Although the original version of DistMesh required actual signed distance functions for the geometry, this condition was relaxed in later versions. The algorithm actually accepts any smooth function with an implicitly defined zero level-set.

This can be very convenient when defining non-trivial geometries. For example, an ellipse with semi-axes 2 and 1 can now be represented as the zero level-set of $d(x,y)=(x/2)^2 + (y/1)^2 - 1$. Note that this is not the actual Euclidean distance function for the ellipse.

```@example introduction
fd(p) = (p[1]/2)^2 + (p[2]/1)^2 - 1
bbox = ((-2,-1), (2,1))
msh = distmesh2d(fd, huniform, 0.2, bbox)
plot(msh)

```

### Non-uniform size function

For non-uniform element sizes, we provide a (relative) size function $h(x,y)$. In general, it is best to make this an absolute size function that gives the actual desired edge lengths at point . To achieve this, you set the initial edge lengths `hmin` to the smallest value of $h(x,y)$ in the domain.

For example, to set the element sizes to  at a source and let the element sizes increase linearly away from the source, we use a size function $h(x,y)=h_\mathrm{min} + 0.3d(x,y)$ where $d(x,y)$ is the distance function of the source.

```@example introduction
# Point source
hmin = 0.01
fh(p) = hmin + 0.3 * dcircle(p, r=0)
msh = distmesh2d(dcircle, fh, hmin, ((-1,-1), (1,1)))
plot(msh)

```

Multiple size function sources can be combined using the `min` function:

```@example introduction
# Point and line sources
fd(p) = drectangle(p, 0, 1, 0, 1)
fh(p) = min(min(0.01 + 0.3*abs(dcircle(p, r=0)),
                0.025 + 0.3*abs(dpoly(p, [(0.3,0.7), (0.7,0.5)]))),
            0.15)
# Note: we explicitly add corners to pfix for the rectangle
pfix = ((0,0), (1,0), (0,1), (1,1))
msh = distmesh2d(fd, fh, 0.01, ((0,0), (1,1)), pfix)
plot(msh)

```

### Randomness and Reproducibility

The DistMesh algorithm is notoriously *chaotic*, or sensitive to initial conditions. This means that very small perturbations in the mesh calculations can grow to large changes in the mesh (resulting in completely different meshes).

However, running on the same deterministic computer, if you repeat a `distmesh2d` call twice with uniform size functions, it usually gives two identical meshes because it uses a deterministic grid-based initialization:

```@example introduction
msh1 = distmesh2d(dcircle, huniform, 0.2, ((-1,-1), (1,1)))
msh2 = distmesh2d(dcircle, huniform, 0.2, ((-1,-1), (1,1)))
msh1.p == msh2.p && msh1.t == msh2.t

```

For **non-uniform size functions**, DistMesh uses the `rand` function for the initial point distribution (rejection sampling). This means two meshes with identical inputs will in general not be the same:

```@example introduction
fh(p) = 0.01 + 0.3 * dcircle(p, r=0)
msh1 = distmesh2d(dcircle, fh, 0.01, ((-1,-1), (1,1)))

```

```@example introduction
msh2 = distmesh2d(dcircle, fh, 0.01, ((-1,-1), (1,1)))

```

```@example introduction
# Check if they are identical (likely false)
msh1.p == msh2.p && msh1.t == msh2.t

```

Many times this is undesirable (e.g., for regression testing or debugging). You can seed the random number generator to make the meshes identical:

```@example introduction
using Random
fh(p) = 0.01 + 0.3 * dcircle(p, r=0)

Random.seed!(1234)
msh1 = distmesh2d(dcircle, fh, 0.01, ((-1,-1), (1,1)))

```

Now we generate the second mesh with the same seed:

```@example introduction
Random.seed!(1234)
msh2 = distmesh2d(dcircle, fh, 0.01, ((-1,-1), (1,1)))

```

Finally, we verify they are identical:

```@example introduction
msh1.p == msh2.p && msh1.t == msh2.t

```

### Save/export meshes

DistMesh does not provide built-in mesh export functionality to external formats. However, it is easy to write specialized routines, e.g., based on text files.

A common format is to have an initial line with metadata (number of nodes and elements), followed by comma-separated rows for node positions and element connectivities:

```@example introduction
using DelimitedFiles

function savemesh(msh::DMesh, fname)
    open(fname, "w") do io
        p, t = as_arrays(msh)
        println(io, "$(length(p)) $(length(t)) # nbr_nodes nbr_elems")
        writedlm(io, p', ',')
        writedlm(io, t', ',')
    end
end

msh = distmesh2d(dcircle, huniform, 0.2, ((-1,-1), (1,1)))
fname = joinpath(tempdir(), "mesh.dat")
savemesh(msh, fname)
println("Mesh saved to file $fname")

```

---

## More about Distance Functions

TODO

---

## More about Size Functions

TODO
