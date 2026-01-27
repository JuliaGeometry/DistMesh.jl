# -------------------------------------------------------------------------
# Mesh Data Structures and Helpers
# -------------------------------------------------------------------------

"""
    DMesh{D, T, N, I}

A lightweight container for mesh data.
- `D`: Spatial dimension (e.g., 2 for 2D coordinates).
- `T`: Floating point type for coordinates (e.g., Float64).
- `N`: Number of vertices per element (e.g., 3 for triangles).
- `I`: Integer type for indices (e.g., Int, Int32).
"""
struct DMesh{D, T, N, I <: Integer}
    p::Vector{SVector{D, T}}
    t::Vector{SVector{N, I}}
end

# Make DMesh iterable so it acts like (p, t)
Base.iterate(m::DMesh, state=1) = iterate((m.p, m.t), state)
Base.eltype(::Type{DMesh{D,T,N,I}}) where {D,T,N,I} = Union{Vector{SVector{D,T}}, Vector{SVector{N,I}}}
Base.length(::DMesh) = 2

"""
    p_view, t_view = as_arrays(m::DMesh)

Return zero-copy views of the mesh nodes and elements.

Note: The shape is `(D x NumPoints)` and `(N x NumElements)`. 
This corresponds to Julia's column-major memory layout (columns are points).
Modifying these arrays will modify the underlying `DMesh`.
"""
function as_arrays(m::DMesh{D,T,N,I}) where {D,T,N,I}
    p_view = reinterpret(reshape, T, m.p)
    t_view = reinterpret(reshape, I, m.t)
    return p_view, t_view
end


# --- Helper for element names ---
element_name(::Val{3}) = "triangle"
element_name(::Val{4}) = "tetrahedron"
element_name(::Val{N}) where N = "$N-simplex"

# 1. Compact Show (Standard)
function Base.show(io::IO, m::DMesh{D,T,N,I}) where {D,T,N,I}
    print(io, "DMesh{$D,$T}($(length(m.p))n, $(length(m.t))e)")
end

# 2. Rich Show (MIME)
function Base.show(io::IO, ::MIME"text/plain", m::DMesh{D,T,N,I}) where {D,T,N,I}
    print(io, "DMesh: $(D)D, ")
    print(io, "$(length(m.p)) nodes, ")
    print(io, "$(length(m.t)) $(element_name(Val(N))) elements")
end


# Global flag to track if we have warned the user yet
const _has_warned_plot = Ref(false)

function live_plot(args...)
    if !_has_warned_plot[]
        @warn "Live plotting was requested, but no plotting backend is loaded. Try `using Plots`."
        _has_warned_plot[] = true
    end
    return nothing
end
