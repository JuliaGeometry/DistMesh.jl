# -------------------------------------------------------------------------
# Mesh Data Structures and Helpers
# -------------------------------------------------------------------------

"""
    DMesh{D, T, G, N, I}

A lightweight container for mesh data.
- `D`: Spatial dimension (e.g., 2 for 2D coordinates).
- `T`: Floating point type for coordinates (e.g., Float64).
- `G`: Element geometry type (e.g., Simplex{2} or Block{3}).
- `N`: Number of vertices per element (e.g., 3 for triangles).
- `I`: Integer type for indices (e.g., Int, Int32).
"""
struct DMesh{D, T, G <: ElementGeometry, N, I <: Integer}
    p::Vector{SVector{D, T}}
    t::Vector{SVector{N, I}}

    # Inner constructor strictly enforces N matches the geometry
    function DMesh(p::Vector{SVector{D, T}}, t::Vector{SVector{N, I}}, geom::G) where {D, T, N, I, G <: ElementGeometry}
        @assert N == nvertices(geom) "Mismatch: ElementGeometry expects $(nvertices(geom)) nodes, but elements have $N nodes."
        new{D, T, G, N, I}(p, t)
    end
end

# -------------------------------------------------------------------------
# Constructors
# -------------------------------------------------------------------------

# 1. Outer constructor: Deduce geometry from SVector inputs
function DMesh(p::Vector{SVector{D, T}}, t::Vector{SVector{N, I}}) where {D, T, N, I <: Integer}
    return DMesh(p, t, find_elgeom(D, N))
end

# 2. Outer constructor: Explicit geometry with NTuple inputs
function DMesh(p::AbstractVector{<:NTuple{D, T}}, t::AbstractVector{<:NTuple{N, I}}, geom::ElementGeometry) where {D, T, N, I <: Integer}
    return DMesh(SVector{D, T}.(p), SVector{N, I}.(t), geom)
end

# 3. Outer constructor: Deduce geometry from NTuple inputs
function DMesh(p::AbstractVector{<:NTuple{D, T}}, t::AbstractVector{<:NTuple{N, I}}) where {D, T, N, I <: Integer}
    return DMesh(SVector{D, T}.(p), SVector{N, I}.(t), find_elgeom(D, N))
end

# -------------------------------------------------------------------------
# Display
# -------------------------------------------------------------------------

# Hide the type complexity in the REPL
function Base.show(io::IO, m::DMesh{D, T, G, N, I}) where {D, T, G, N, I}
    geom_name = name(G())
    print(io, "$(D)D DMesh ($geom_name) with $(length(m.p)) points and $(length(m.t)) elements")
end

# -------------------------------------------------------------------------
# Convertors
# -------------------------------------------------------------------------

# Make DMesh iterable so it acts like (p, t)
Base.iterate(m::DMesh, state=1) = iterate((m.p, m.t), state)
Base.eltype(::Type{DMesh{D, T, G, N, I}}) where {D, T, G, N, I} = Union{Vector{SVector{D, T}}, Vector{SVector{N, I}}}
Base.length(::DMesh) = 2

"""
    p_view, t_view = as_arrays(m::DMesh)

Return zero-copy views of the mesh nodes and elements.

The shape is `(D x NumPoints)` and `(N x NumElements)`. 
Modifying these arrays will modify the underlying `DMesh`.
"""
function as_arrays(m::DMesh{D, T, G, N, I}) where {D, T, G, N, I}
    p_view = reinterpret(reshape, T, m.p)
    t_view = reinterpret(reshape, I, m.t)
    return p_view, t_view
end

# -------------------------------------------------------------------------
# Other utilities
# -------------------------------------------------------------------------

# Global flag to track if we have warned the user yet
const _has_warned_plot = Ref(false)

function live_plot(args...)
    if !_has_warned_plot[]
        @warn "Live plotting was requested, but no plotting backend is loaded. Try `using Plots`."
        _has_warned_plot[] = true
    end
    return nothing
end
