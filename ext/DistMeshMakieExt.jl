module DistMeshMakieExt

using DistMesh
import Makie

const MESH_COLOR = "#DDEEFF"

function Makie.plot(m::DMesh{2}; args...)
    f = Makie.Figure()
    ax = Makie.Axis(f[1,1], aspect=Makie.DataAspect())

    # Build a vector of Polygons — works for any element shape (triangles, quads, ...)
    polys = [Makie.Polygon([Makie.Point2f(m.p[i]) for i in el]) for el in m.t]
    Makie.poly!(ax, polys, color=MESH_COLOR, strokewidth=1)
    return f
end

const _has_warned_live = Ref(false)

function DistMesh.live_plot(m::DMesh)
    if !_has_warned_live[]
        @warn "Live plotting is not supported for this Makie backend.\n" *
              "Switch to `using GLMakie` for animations."
        _has_warned_live[] = true
    end
    return nothing
end

end
