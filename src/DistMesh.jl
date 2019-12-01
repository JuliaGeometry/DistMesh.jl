module DistMesh

using GeometryBasics,
      LinearAlgebra,
      TetGen

abstract type AbstractDistMeshAlgorithm end

"""
    DistMeshSetup

    The original DistMesh retriangulation/termination criterea. Will trigger a
    retriangulation once a point movement greater than ttol is exceeded.
    Will terminate once an iteration produces a point movement less than ptol.

    Takes Keyword arguments as follows:

    iso (default: 0): Value of which to extract the isosurface, inside surface is negative
    deltat (default: 0.1): the fraction of edge displacement to apply each iteration
    sort (default:false): If true and no fixed points, sort points using a hilbert sort.
    sort_interval (default:20) Use hilbert sort after the specified retriangulations
    distribution (default: :regular) Initial point distribution, either :regular or :packed.
"""
struct DistMeshSetup{T} <: AbstractDistMeshAlgorithm
    iso::T
    deltat::T
    ttol::T
    ptol::T
    sort::Bool # use hilbert sort to cache-localize points
    sort_interval::Int # retriangulations before resorting
    nonlinear::Bool # uses nonlinear edge force
    distribution::Symbol # intial point distribution
end

function DistMeshSetup(;iso=0,
                        ptol=.001,
                        deltat=0.05,
                        ttol=0.02,
                        sort=false,
                        sort_interval=20,
                        nonlinear=false,
                        distribution=:regular)
    T = promote_type(typeof(iso),typeof(ptol),typeof(deltat), typeof(ttol))
    DistMeshSetup{T}(iso,
                     deltat,
                     ttol,
                     ptol,
                     sort,
                     sort_interval,
                     nonlinear,
                     distribution)
end

"""
    DistMeshQuality

    iso (default: 0): Value of which to extract the iso surface, inside negative
    deltat (default: 0.1): the fraction of edge displacement to apply each iteration
"""
struct DistMeshQuality{T} <: AbstractDistMeshAlgorithm
    iso::T
    deltat::T
    minimum::T
    sort::Bool # use hilbert sort to cache-localize points
    sort_interval::Int # retriangulations before resorting
    nonlinear::Bool # uses nonlinear edge force
    distribution::Symbol # intial point distribution
end

function DistMeshQuality(;iso=0,
                        deltat=0.05,
                        minimum=0.65,
                        sort=false,
                        sort_interval=20,
                        nonlinear=false,
                        distribution=:regular)
    T = promote_type(typeof(iso),typeof(deltat))
    DistMeshQuality{T}(iso,
                     deltat,
                     minimum,
                     sort,
                     sort_interval,
                     nonlinear,
                     distribution)
end

"""
    DistMeshStatistics

        Statistics about the convergence between iterations
"""
struct DistMeshStatistics{T}
    maxmove::Vector{T} # max point move in an iteration
    maxdp::Vector{T} # max displacmeent induced by an edge
    average_qual::Vector{T}
    median_qual::Vector{T}
    minimum_qual::Vector{T}
    maximum_qual::Vector{T}
    retriangulations::Vector{Int} # Iteration num where retriangulation occured
end

DistMeshStatistics() = DistMeshStatistics{Float64}([],[],[],[],[],[],[])

"""
Uniform edge length specifier.
"""
struct HUniform end

include("common.jl")
include("diff.jl")
include("pointdistribution.jl")
include("distmeshnd.jl")
include("tetgen.jl")
include("quality.jl")
include("decompositions.jl")
include("hilbertsort.jl")

#export distmeshsurface
export distmesh, DistMeshSetup, DistMeshStatistics, HUniform

end # module
