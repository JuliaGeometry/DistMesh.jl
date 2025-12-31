# ## Polygon Mesh
#
# This example demonstrates how to mesh a custom polygon using `dpoly`.

using DistMesh
using Plots

# ### 1. Define Vertices
pv = [(-0.4, -0.5), (0.4, -0.2), (0.4, -0.7), (1.5, -0.4),
      (0.9, 0.1), (1.6, 0.8), (0.5, 0.5), (0.2, 1.0),
      (0.1, 0.4), (-0.7, 0.7), (-0.4, -0.5)]

# ### 2. Create Distance Function
# `dpoly` computes the distance from points `p` to the polygon defined by `pv`.
# We wrap it in a new function because `distmesh2d` expects a function of `p`.
fd(p) = dpoly(p, pv)

# ### 3. Generate Mesh
# We use a uniform mesh size of h0 = 0.1. 
# The bounding box must enclose the polygon.

bbox = ((-1,-1), (2,1))
h0 = 0.15

msh = distmesh2d(fd, huniform, h0, bbox, pv)

# ### 4. Visualize
# Plot the resulting mesh.

plot(msh)
