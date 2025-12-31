module DistMeshGLMakieExt

using DistMesh
using GLMakie

const meshgreen = Plots.RGBX(0.8, 1, 0.8)

function GLMakie.plot(m::DMesh{2,T,N,I};
    labels=(), reltol=1e-3, abstol=Inf, maxref=6,
    colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
    # colors: elements, int edges, bnd edges, ho-nodes, vertices
) where {T,N,I}

    # elem_lines, int_lines, bnd_lines, elem_mid =
    #     viz_mesh(m, reltol=reltol, abstol=abstol, maxref=maxref)

    # h = Plots.plot(aspect_ratio=:equal, leg=false)
    # Plots.plot!(elem_lines[:, 1], elem_lines[:, 2], seriestype=:shape,
    #             linecolor=nothing, fillcolor=colors[1])
    # Plots.plot!(int_lines[:, 1], int_lines[:, 2], linewidth=1, color=colors[2])
    # Plots.plot!(bnd_lines[:, 1], bnd_lines[:, 2], linewidth=2, color=colors[3])

    # if isa(labels, Symbol)
    #     labels = (labels,)
    # end

    # if :nodes in labels
    #     scatter!(m.x[:,1], m.x[:,2], marker=(:circle, 2, 1.0, colors[4]))
    # end

    # if :elements in labels
    #     for (iel,mid) in enumerate(eachrow(elem_mid))
    #         annotate!(mid[1], mid[2], text("$iel", 8, :center))
    #     end
    # end

    # return h
end

end
