module DistMeshCairoMakieExt

using DistMesh
using CairoMakie

const MESH_COLOR = "#DDEEFF"

function CairoMakie.plot(m::DMesh{2}; args...) 
    f = Figure(size=(800, 800))
    ax = Axis(f[1,1], aspect=DataAspect())
    
    p, t = as_arrays(m)
    poly!(ax, p', t', color=MESH_COLOR, strokewidth=1)
    return f
end

const _has_warned_live = Ref(false)

function DistMesh.live_plot(m::DMesh)
    if !_has_warned_live[]
        @warn "Live plotting disabled for CairoMakie. \n" *
              "CairoMakie is for static plots (PDF/PNG). \n" *
              "Switch to `using GLMakie` for animations."
        _has_warned_live[] = true
    end
    return nothing
end

end

