########################################################################
# Internal Type Aliases

const Point2d = SVector{2, Float64}
const Point3d = SVector{3, Float64}
const Index2 = SVector{2, Int32}   # For Edges
const Index3 = SVector{3, Int32}   # For Triangles

# -------------------------------------------------------------------------
# Mesh Data Structures and Helpers
# -------------------------------------------------------------------------

"""
    DMesh{D, T, G, N, I}

A lightweight container for mesh data.
- `D`: Spatial dimension (e.g., 2 for 2D coordinates).
- `T`: Floating point type for coordinates (e.g., Float64).
- `G`: Element topology type (e.g., Simplex{2} or Block{3}).
- `N`: Number of vertices per element (e.g., 3 for triangles).
- `I`: Integer type for indices (e.g., Int, Int32).
"""
struct DMesh{D, T, G <: ElementTopology, N, I <: Integer}
    p::Vector{SVector{D, T}}
    t::Vector{SVector{N, I}}

    # Inner constructor strictly enforces N matches the topology
    function DMesh(p::Vector{SVector{D, T}}, t::Vector{SVector{N, I}}, geom::G) where {D, T, N, I, G <: ElementTopology}
        @assert N == nvertices(geom) "Mismatch: ElementTopology expects $(nvertices(geom)) nodes, but elements have $N nodes."
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
function DMesh(p::AbstractVector{<:NTuple{D, T}}, t::AbstractVector{<:NTuple{N, I}}, geom::ElementTopology) where {D, T, N, I <: Integer}
    return DMesh(SVector{D, T}.(p), SVector{N, I}.(t), geom)
end

# 3. Outer constructor: Deduce geometry from NTuple inputs
function DMesh(p::AbstractVector{<:NTuple{D, T}}, t::AbstractVector{<:NTuple{N, I}}) where {D, T, N, I <: Integer}
    return DMesh(SVector{D, T}.(p), SVector{N, I}.(t), find_elgeom(D, N))
end

# -------------------------------------------------------------------------
# Matrix Constructors
# -------------------------------------------------------------------------

# Helper function to bridge runtime matrix dimensions to compile-time SVector parameters
function _mat_to_svec(mat::AbstractMatrix{T}, ::Val{K}) where {T, K}
    dense_mat = mat isa Matrix ? mat : Matrix(mat)
    return copy(reinterpret(reshape, SVector{K, T}, dense_mat))
end

# 4. Outer constructor: Matrices D-by-NP and N-by-NT
function DMesh(p::AbstractMatrix{T}, t::AbstractMatrix{I}, geom::ElementTopology) where {T, I <: Integer}
    D = size(p, 1)
    N = size(t, 1)
    
    p_vec = _mat_to_svec(p, Val(D))
    t_vec = _mat_to_svec(t, Val(N))
    
    return DMesh(p_vec, t_vec, geom)
end

function DMesh(p::AbstractMatrix{T}, t::AbstractMatrix{I}) where {T, I <: Integer}
    D = size(p, 1)
    N = size(t, 1)
    
    return DMesh(p, t, find_elgeom(D, N))
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
