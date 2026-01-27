module DistMeshPlotsExt

using DistMesh
using Plots

const MESH_COLOR = Plots.RGBX(0.8, 0.9, 1.0)

function Plots.plot(m::DMesh{2,T,N,I}; args...) where {T,N,I}
    pxy(d) = [ tix==0 ? NaN : m.p[tt[tix]][d] for tix in (1,2,3,1,0), tt in m.t ]
    args = (args...,
            seriestype=:shape,
            aspect_ratio=:equal,
            leg=false,
            fillcolor=MESH_COLOR)
    h = Plots.plot(vec(pxy(1)), vec(pxy(2)); args...)
end

DistMesh.live_plot(m::DMesh; args...) = display(plot(m))

end
