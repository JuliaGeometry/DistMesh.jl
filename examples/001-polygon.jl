# ## Polygon Mesh

using DistMesh
using GLMakie

# First, define the vertices of the polygon as a list of points (tuples)
pv = [(-0.4, -0.5), (0.4, -0.2), (0.4, -0.7), (1.5, -0.4),
      (0.9, 0.1), (1.6, 0.8), (0.5, 0.5), (0.2, 1.0),
      (0.1, 0.4), (-0.7, 0.7), (-0.4, -0.5)];

# Next, use the `dpoly` function and these vertices to define the distance function
fd(p) = dpoly(p, pv)

# We use a uniform mesh size function, with a given desired edge length
hmin = 0.15
fh = huniform

# DistMesh also needs a bounding box, which must enclose the polygon
bbox = ((-1,-1), (2,1))

# Sharp corners have to be provided in the `pfix` argument
pfix = pv

# Generate the mesh
msh = distmesh2d(fd, fh, hmin, bbox, pfix)

# Finally, plot the resulting mesh
plot(msh)
