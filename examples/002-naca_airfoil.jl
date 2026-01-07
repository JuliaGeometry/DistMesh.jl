# ## NACA Airfoil Mesh

using DistMesh
using GLMakie
using StaticArrays

const Point2d = SVector{2, Float64}

function hnaca(p; hlead=0.01, htrail=0.04, hmax=2.0)
    minimum((hlead + 0.3 * dcircle(p, c=(0, 0), r=0),
        htrail + 0.3 * dcircle(p, c=(1, 0), r=0),
        hmax))
end

function fixnaca(; htrail=0.04)
    a0, a14... = naca_coeffs
    fixx = 1 .- htrail * cumsum(1.3 .^ (0:4))
    fixy = a0 * sqrt.(fixx) .+ fixx .^ (1:4)' * a14
    fix = vcat(Point2d[(0,0),(1,0)], [ Point2d[(x,y),(x,-y)] for (x,y) in zip(fixx,fixy) ]...)
end

function dm_naca(; hlead=0.01, htrail=0.04, hmax=2.0, circx=2.0, circr=4.0)
    dfcn(p) = ddiff(dcircle(p, c=(circx, 0), r=circr), dnaca(p))
    hfcn(p) = hnaca(p; hlead=hlead, htrail=htrail, hmax=hmax)

    fix = Point2d[(1,0),(0,1),(-1,0),(0,-1)] .* circr .+ Point2d[(circx,0)]
    fix = vcat(fix, fixnaca(htrail=htrail))
    bbox = [(circx - circr, -circr), (circx + circr, circr)]
    h0 = minimum((hlead, htrail, hmax))
    dfcn, hfcn, h0, bbox, fix
end

msh = distmesh2d(dm_naca()...)

plot(msh)
