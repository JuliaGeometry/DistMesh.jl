function live_plot(args...)
    @warn "Live plotting was requested, but no plotting backend is loaded. Try `using GLMakie`." maxlog=1
    return nothing
end

"""
    get_camera_view(ax; digits=1)

Extract and print the current camera position of a Makie LScene.
Note: **Requires Makie (or GLMakie/WGLMakie) to be loaded.**
"""
function get_camera_view(args...; kwargs...)
    error("`get_camera_view` requires Makie. Please run `using GLMakie` (or Makie) first.")
end
