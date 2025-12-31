include("distmesh.jl")

function ex5polygon(; h0=0.1)
    pv = SF2[(-0.4, -0.5), (0.4, -0.2), (0.4, -0.7), (1.5, -0.4),
             (0.9, 0.1), (1.6, 0.8), (0.5, 0.5), (0.2, 1.0),
             (0.1, 0.4), (-0.7, 0.7), (-0.4, -0.5)]
    dfcn(p) = dpoly(p, pv)
    dfcn, huniform, h0, [(-1, -1), (2, 1)], pv
end

function hnaca(p; hlead=0.01, htrail=0.04, hmax=2.0)
    minimum((hlead + 0.3 * dcircle(p, c=(0, 0), r=0),
        htrail + 0.3 * dcircle(p, c=(1, 0), r=0),
        hmax))
end

function fixnaca(; htrail=0.04)
    a0, a14... = naca_coeffs
    fixx = 1 .- htrail * cumsum(1.3 .^ (0:4))
    fixy = a0 * sqrt.(fixx) .+ fixx .^ (1:4)' * a14
    fix = vcat(SF2[(0,0),(1,0)], [ SF2[(x,y),(x,-y)] for (x,y) in zip(fixx,fixy) ]...)
end

function ex6naca(; hlead=0.01, htrail=0.04, hmax=2.0, circx=2.0, circr=4.0)
    dfcn(p) = ddiff(dcircle(p, c=(circx, 0), r=circr), dnaca(p))
    hfcn(p) = hnaca(p; hlead=hlead, htrail=htrail, hmax=hmax)

    fix = SF2[(1,0),(0,1),(-1,0),(0,-1)] .* circr .+ SF2[(circx,0)]
    fix = vcat(fix, fixnaca(htrail=htrail))
    bbox = [(circx - circr, -circr), (circx + circr, circr)]
    h0 = minimum((hlead, htrail, hmax))
    dfcn, hfcn, h0, bbox, fix
end

for ex in (ex5polygon, ex6naca)
    p, t = distmesh2d(ex()...)
    display(simpplot(p, t))
end

using Random
#for s = 1:1000
#    @show s
#    Random.seed!(s)
#    p,t = distmesh2d(ex6naca()...);
#    @show minimum(simpqual(p,t))
#end

#Random.seed!(9)
#distmesh2d(ex6naca()...)
