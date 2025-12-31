module DistMeshPlotsExt

using DistMesh
using Plots

const meshgreen = Plots.RGBX(0.8, 1, 0.8)

function Plots.plot(m::DMesh{2,T,N,I}; args...) where {T,N,I}
    pxy(d) = [ tix==0 ? NaN : m.p[tt[tix]][d] for tix in (1,2,3,1,0), tt in m.t ]
    args = (args...,
            seriestype=:shape,
            aspect_ratio=:equal,
            leg=false,
            fillcolor=Plots.RGB(0.8, 0.9, 1.0))
    h = Plots.plot(vec(pxy(1)), vec(pxy(2)); args...)
end

end
