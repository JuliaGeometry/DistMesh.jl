module DistMeshMakieExt

using DistMesh
using Makie
import GeometryBasics

const MESH_COLOR = "#DDEEFF"

# ---------------------------------------------------------
# 1. Helpers
# ---------------------------------------------------------

"""
    get_canvas(; aspect=DataAspect())

Gets the current figure/axis if they exist, or creates new ones.
"""
function get_canvas(; aspect=DataAspect())
    fig = current_figure()
    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1], aspect=aspect)
        return fig, ax
    else
        ax = current_axis()
        if isnothing(ax)
            ax = Axis(fig[1,1], aspect=aspect)
        end
        return fig, ax
    end
end

"""
    to_makie_mesh(m::DMesh)

Converts a DMesh into a GeometryBasics.normal_mesh.
This strict conversion avoids internal Makie conversions, 
yielding maximum performance for updates across ALL backends.
"""
function to_makie_mesh(m::DMesh{2,T,Simplex{2}}) where {T}
    p, t = as_arrays(m)
    pts = GeometryBasics.Point2f[GeometryBasics.Point2f(col) for col in eachcol(p)]
    faces = GeometryBasics.GLTriangleFace[GeometryBasics.GLTriangleFace(col...) for col in eachcol(t)]
    return GeometryBasics.normal_mesh(GeometryBasics.Mesh(pts, faces))
end

function to_makie_mesh(m::DMesh{3,T,Simplex{2}}) where {T}
    pts = GeometryBasics.Point3f.(m.p)
    faces = GeometryBasics.GLTriangleFace[GeometryBasics.GLTriangleFace(f...) for f in m.t]
    return GeometryBasics.normal_mesh(GeometryBasics.Mesh(pts, faces))
end

# Quad mesh: split each quad (CCW: 1-2-3-4) into two triangles (1,2,3) and (1,3,4)
function to_makie_mesh(m::DMesh{2,T,Block{2}}) where {T}
    p, t = as_arrays(m)
    pts = GeometryBasics.Point2f[GeometryBasics.Point2f(col) for col in eachcol(p)]
    faces = GeometryBasics.GLTriangleFace[]
    sizehint!(faces, 2 * size(t, 2))
    for col in eachcol(t)
        a, b, c, d = col
        push!(faces, GeometryBasics.GLTriangleFace(a, b, c))
        push!(faces, GeometryBasics.GLTriangleFace(a, c, d))
    end
    return GeometryBasics.normal_mesh(GeometryBasics.Mesh(pts, faces))
end

# ---------------------------------------------------------
# 2. Standard Plot (Static)
# ---------------------------------------------------------

# Fast path for 2D triangle meshes
function Makie.plot(m::DMesh{2,T,Simplex{2}}; args...) where {T}
    f, ax = get_canvas()
    empty!(ax)
    fast_mesh = to_makie_mesh(m)
    poly!(ax, fast_mesh, color=MESH_COLOR, strokewidth=1)
    return f
end

# Generic fallback for 2D polygon meshes (Quads, etc.)
function Makie.plot(m::DMesh{2}; args...)
    f, ax = get_canvas()
    empty!(ax)
    polys = [GeometryBasics.Polygon([GeometryBasics.Point2f(m.p[i]) for i in el]) for el in m.t]
    poly!(ax, polys, color=MESH_COLOR, strokewidth=1)
    return f
end

# ---------------------------------------------------------
# 3. Live Plot (Dynamic / Animation)
# ---------------------------------------------------------

# Returns true if the mesh's bounding box has moved outside the current view,
# or the view has grown more than 2x larger than the mesh in any dimension.
# Cheap enough to run every frame.
function _needs_relimit(m::DMesh{2}, ax)
    p, _ = as_arrays(m)
    xmin, xmax = extrema(view(p, 1, :))
    ymin, ymax = extrema(view(p, 2, :))
    lims = ax.finallimits[]
    lx, ly = lims.origin[1], lims.origin[2]
    rx, ry = lx + lims.widths[1], ly + lims.widths[2]
    return xmin < lx || xmax > rx || ymin < ly || ymax > ry ||
           (rx - lx) > 2 * (xmax - xmin) || (ry - ly) > 2 * (ymax - ymin)
end

function DistMesh.live_plot(m::DMesh{2,T,Simplex{2}}) where {T}
    # Check if the active backend is interactive
    backend_name = string(Makie.current_backend())
    if !occursin("GLMakie", backend_name) && !occursin("WGLMakie", backend_name)
        @warn "Live plotting requires an interactive backend. Switch to `using GLMakie` for animations." maxlog=1
    end

    f, ax = get_canvas()
    fast_mesh = to_makie_mesh(m)

    if isempty(ax.scene.plots)
        poly!(ax, fast_mesh, color=MESH_COLOR, strokewidth=1)
        autolimits!(ax)
        display(f)
    else
        plt = ax.scene.plots[1]
        plt[1][] = fast_mesh
        _needs_relimit(m, ax) && autolimits!(ax)
    end

    # sleep(0.01) # You might even want to skip the sleep if it's not GLMakie to speed up the dummy loop!
    if occursin("GLMakie", backend_name) || occursin("WGLMakie", backend_name)
        sleep(0.01)
    end

    return f
end

function DistMesh.live_plot(m::DMesh{2,T,Block{2}}) where {T}
    backend_name = string(Makie.current_backend())
    if !occursin("GLMakie", backend_name) && !occursin("WGLMakie", backend_name)
        @warn "Live plotting requires an interactive backend. Switch to `using GLMakie` for animations." maxlog=1
    end

    f, ax = get_canvas()
    polys = [GeometryBasics.Polygon([GeometryBasics.Point2f(m.p[i]) for i in el]) for el in m.t]

    if isempty(ax.scene.plots)
        poly!(ax, polys, color=MESH_COLOR, strokewidth=1)
        autolimits!(ax)
        display(f)
    else
        # poly!'s input observable is typed as Vector{Polygon}, so in-place update works.
        plt = ax.scene.plots[1]
        plt[1][] = polys
        _needs_relimit(m, ax) && autolimits!(ax)
    end

    if occursin("GLMakie", backend_name) || occursin("WGLMakie", backend_name)
        sleep(0.01)
    end

    return f
end

# ---------------------------------------------------------
# 4. 3D Utilities & Plotting
# ---------------------------------------------------------

function DistMesh.get_camera_view(ax::LScene; digits=1)
    cam = Makie.cameracontrols(ax.scene)
    eye = round.(Tuple(cam.eyeposition[]), digits = digits)
    look = round.(Tuple(cam.lookat[]), digits = digits)
    up = round.(Tuple(cam.upvector[]), digits = digits)
    println("_campos = (Vec3f$eye, Vec3f$look, Vec3f$up)")
    println("update_cam!(ax.scene, _campos...)")
end

function Makie.plot(m::DMesh{3,T,Simplex{2}}; 
                      color=("#CCFFCC", "#FFAAAA", "#DDEEFF"),
                      campos=nothing, elems=nothing, 
                      plot_rest=true, rest_color=(:lightgray, 0.1)) where {T}

    fig = current_figure()
    local ax
    local existing_campos = nothing
    
    if isnothing(fig) || isempty(fig.content)
        fig = Figure()
        ax = LScene(fig[1, 1]) 
    else
        ax = fig.content[1]
        cam = Makie.cameracontrols(ax.scene)
        existing_campos = (cam.eyeposition[], cam.lookat[], cam.upvector[])
        
        for p in copy(ax.scene.plots)
            if p isa Makie.Mesh || p isa Makie.Wireframe
                delete!(ax.scene, p)
            end
        end
    end

    fast_mesh = to_makie_mesh(m)
    pts = GeometryBasics.coordinates(fast_mesh)
    fs = GeometryBasics.faces(fast_mesh)

    elems === nothing && (elems = 1:length(fs))
    elems isa Tuple || (elems = (elems,))
    color isa Tuple || (color = (color,))
    length(color) == 1 && (color = (color...,))
    selected_indices = union(elems...) 

    if plot_rest
        rest_indices = setdiff(1:length(fs), selected_indices)
        if !isempty(rest_indices)
            rest_m = GeometryBasics.Mesh(pts, fs[rest_indices])
            mesh!(ax, rest_m, color=rest_color, transparency=true)
            wireframe!(ax, rest_m, color=(:black, 0.1), linewidth=1.0, transparency=true)
        end
    end

    for (idxs, c) in zip(elems, color)
        sub_faces = fs[idxs]
        sub_m = GeometryBasics.Mesh(pts, sub_faces)
        mesh!(ax, sub_m, color=c)
        wireframe!(ax, sub_m, color=:black, linewidth=1.0)
    end
    
    display(fig)
    
    if campos !== nothing
        update_cam!(ax.scene, campos...)
    elseif existing_campos !== nothing
        update_cam!(ax.scene, existing_campos...)
    end
    
    return fig, ax
end

end
