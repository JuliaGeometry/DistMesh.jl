module DistMeshGLMakieExt

using DistMesh
using GLMakie
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
        fig = Figure() # size=(600, 600)
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
    to_gl_mesh(m::DMesh)

Converts a DistMesh DMesh object into a GeometryBasics.normal_mesh
suitable for high-performance GLMakie plotting and updates.
"""
function to_gl_mesh(m::DMesh{2,T,Simplex{2}}) where {T}
    p, t = as_arrays(m)

    # 1. strict conversion to Float32 Points (GL standard)
    pts = GeometryBasics.Point2f[GeometryBasics.Point2f(col) for col in eachcol(p)]
    
    # 2. strict conversion to standard Triangle Faces
    faces = GeometryBasics.GLTriangleFace[GeometryBasics.GLTriangleFace(col...) for col in eachcol(t)]
    
    # 3. Create Mesh and compute Normals
    #    (Required because Makie initializes with normals by default; 
    #     we must match that type for updates to work)
    return GeometryBasics.normal_mesh(GeometryBasics.Mesh(pts, faces))
end

# ---------------------------------------------------------
# 2. Standard Plot (Static)
# ---------------------------------------------------------

# Fast GL path for triangle meshes
function GLMakie.plot(m::DMesh{2,T,Simplex{2}}; args...) where {T}
    f, ax = get_canvas()
    empty!(ax)
    gl_mesh = to_gl_mesh(m)
    poly!(ax, gl_mesh, color=MESH_COLOR, strokewidth=1)
    return f
end

# ---------------------------------------------------------
# 3. Live Plot (Dynamic / Animation)
# ---------------------------------------------------------

function DistMesh.live_plot(m::DMesh{2,T,Simplex{2}}) where {T}
    f, ax = get_canvas()
    
    # Always convert the incoming data to the strict GL mesh format
    gl_mesh = to_gl_mesh(m)

    if isempty(ax.scene.plots)
        # --- INITIALIZATION ---
        poly!(ax, gl_mesh, color=MESH_COLOR, strokewidth=1)
        autolimits!(ax)
        display(f) 
    else
        # --- UPDATE ---
        plt = ax.scene.plots[1]
        plt[1][] = gl_mesh
        # autolimits!(ax)   # technically needed, but makes animation non-smooth
    end
    
    sleep(0.01)
    return f
end

# live_plot fallback for non-triangle meshes
function DistMesh.live_plot(m::DMesh{2})
    @warn "Live plotting is only supported for triangle meshes (Simplex{2})."
    return nothing
end

end
