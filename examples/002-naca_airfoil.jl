# ## NACA Airfoil Mesh
#
# This example shows how to generate a mesh around a NACA0012 airfoil.

# First, load the required packages and define a convenient short-hand
# for static 2D points.

using DistMesh
using GLMakie
using StaticArrays

const Point2d = SVector{2, Float64}

# First, we define the size function for the NACA airfoil. It consists of the minimum of several point sources and constants:
# * A point source at the tip of the airfoil, with size `hlead`
# * A point source at the traling edge of the airfoil, with size `htrail`
# * A maximum element size in the entire domain `hmax`

function hnaca(p; hlead=0.01, htrail=0.04, hmax=2.0)
    minimum((hlead + 0.3 * dcircle(p, c=(0, 0), r=0),
        htrail + 0.3 * dcircle(p, c=(1, 0), r=0),
        hmax))
end

# Next, we define the fix points. This could be as simple as just the trailing edge at `(1,0)`. However, the mesh quality is improved by providing several points along the surface, spread out consistently with the size function above. Therefore, this function also needs the edge length `htrail` from before. We also add the symmetry point `(0,0)`.

function fixnaca(; htrail=0.04)
    a0, a14... = naca_coeffs
    fixx = 1 .- htrail * cumsum(1.3 .^ (0:4))
    fixy = a0 * sqrt.(fixx) .+ fixx .^ (1:4)' * a14
    fix = vcat(Point2d[(0,0),(1,0)], [ Point2d[(x,y),(x,-y)] for (x,y) in zip(fixx,fixy) ]...)
end

# Finally, we can define a function for generating the arguments to `distmesh2d` for generating the NACA mesh. The parameters are the sizes from before, as well as a center point `(circx,0)` and a radius for the far-field circle boundary.

function dm_naca(; hlead=0.01, htrail=0.04, hmax=2.0, circx=2.0, circr=4.0)
    ## Distance function: Difference between the outer circle and the NACA airfoil
    dfcn(p) = ddiff(dcircle(p, c=(circx, 0), r=circr), dnaca(p))

    ## Size function: Given by `hnaca` above
    hfcn(p) = hnaca(p; hlead=hlead, htrail=htrail, hmax=hmax)

    ## Fix points: Add symmetry points for the circle plus the once from `fixnaca`
    fix = Point2d[(1,0),(0,1),(-1,0),(0,-1)] .* circr .+ Point2d[(circx,0)]
    fix = vcat(fix, fixnaca(htrail=htrail))

    ## Bounding box: Determined by the far-field circle
    bbox = [(circx - circr, -circr), (circx + circr, circr)]

    ## The `hmin` parameter needs to be the smallest element size in the domain
    hmin = minimum((hlead, htrail, hmax))
    
    dfcn, hfcn, hmin, bbox, fix
end

# Now we can generate and plot the default mesh
msh = distmesh2d(dm_naca()...)
plot(msh)

# We can modify the mesh, e.g. by shrinking the farfield circle
msh = distmesh2d(dm_naca(circx=1, circr=1.5)...)
plot(msh)

# Finally, we can increase the element sizes at the sources
msh = distmesh2d(dm_naca(circx=1, circr=1.5, hlead=0.05, htrail=0.1)...)
plot(msh)
