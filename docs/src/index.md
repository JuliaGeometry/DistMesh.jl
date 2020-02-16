# DistMesh.jl


DistMesh.jl implements simplex refinement on signed distance functions, or anything that
has a sign, distance, and called like a function. The algorithm was first presented
in 2004 by Per-Olof Persson, and was initially a port of the corresponding Matlab Code.

## What is Simplex Refinement?

In layman's terms, a simplex is either a triangle in the 2D case, or a tetrahedra in the 3D case.

When simulating, you other want a few things from a mesh of simplices:
    - Accurate approximation of boundaries and features
    - Adaptive mesh sizes to improve accuracy
    - Near-Regular Simplicies

DistMesh is designed to address the above.

## Algorithm Overview

The basic processes is as follows:



## Comparison to other refinements

DistMesh generally has a very low memory footprint, and can refine without additional
memory allocation. Similarly, since the global state of simplex qualities is accounted for
in each refinement iteration, this leads to very high quality meshes.

Aside from the above, since DistMesh works on signed distance functions it can handle
complex and varied input data that are not in the form of surface meshes (Piecewise Linear Complicies).

## Difference from the MatLab implementation

Given the same parameters, the Julia implementation of DistMesh will generally perform
4-60 times faster than the MatLab implementation. Delaunay Triangulation in MatLab uses
QHull, whereas DistMesh.jl uses TetGen.

## How do I get a Signed Distance Function?

Here are some libraries that turn gridded and level set data into an approximate signed
distance function:

- Interpolations.jl
- AdaptiveDistanceFields.jl

```@index
```

```@autodocs
Modules = [DistMesh]
```
