module DistMeshGLMakieExt

using DistMesh
using GLMakie
using GeometryBasics 

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
        fig = Figure(size=(800, 800))
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
function to_gl_mesh(m::DMesh{2})
    p, t = as_arrays(m)

    # 1. strict conversion to Float32 Points (GL standard)
    pts = Point2f[Point2f(col) for col in eachcol(p)]
    
    # 2. strict conversion to standard Triangle Faces
    faces = GLTriangleFace[GLTriangleFace(col...) for col in eachcol(t)]
    
    # 3. Create Mesh and compute Normals
    #    (Required because Makie initializes with normals by default; 
    #     we must match that type for updates to work)
    return GeometryBasics.normal_mesh(GeometryBasics.Mesh(pts, faces))
end

# ---------------------------------------------------------
# 2. Standard Plot (Static)
# ---------------------------------------------------------

function GLMakie.plot(m::DMesh{2}; args...) 
    f, ax = get_canvas()
    
    # Standard plot: Clear the axis and draw fresh
    empty!(ax) 
    
    gl_mesh = to_gl_mesh(m)
    poly!(ax, gl_mesh, color=MESH_COLOR, strokewidth=1)
    
    return f
end

# ---------------------------------------------------------
# 3. Live Plot (Dynamic / Animation)
# ---------------------------------------------------------

function DistMesh.live_plot(m::DMesh{2})
    f, ax = get_canvas()
    
    # Always convert the incoming data to the strict GL mesh format
    gl_mesh = to_gl_mesh(m)

    if isempty(ax.scene.plots)
        # --- INITIALIZATION ---
        # If axis is empty, draw for the first time
        poly!(ax, gl_mesh, color=MESH_COLOR, strokewidth=1)
        autolimits!(ax)
    else
        # --- UPDATE ---
        # If plot exists, update the Observable in place (Zero Allocation on GPU)
        plt = ax.scene.plots[1]
        
        # We update the first argument of the plot object
        plt[1][] = gl_mesh
    end
    
    # Yield to the event loop so the window redraws
    sleep(0.01)
    
    return f
end

end
