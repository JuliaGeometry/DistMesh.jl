module DistMesh

using LinearAlgebra,
      TetGen
using GeometryBasics
using GeometryBasics: Triangle, Tetrahedron, Mesh, Polytope, Point

const tetpairs = ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
const tettriangles = ((1,2,3),(1,2,4),(2,3,4),(1,3,4))

abstract type AbstractDistMeshAlgorithm end

"""
    DistMeshSetup

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

Use Tetrahedral quality analysis to control the meshing process

    iso (default: 0): Value of which to extract the iso surface, inside negative
    deltat (default: 0.1): the fraction of edge displacement to apply each iteration
"""
struct DistMeshQuality{T} <: AbstractDistMeshAlgorithm
    iso::T
    deltat::T
    filter_less_than::T # Remove tets less than the given quality
    #allow_n_regressions::Int # Might want this
    termination_quality::T # Once achieved, terminate
    sort::Bool # use hilbert sort to cache-localize points
    sort_interval::Int # retriangulations before resorting
    nonlinear::Bool # uses nonlinear edge force
    distribution::Symbol # intial point distribution
end

function DistMeshQuality(;iso=0,
                        deltat=0.05,
                        filter_less_than=0.02,
                        termination_quality=0.3,
                        sort=false,
                        sort_interval=20,
                        nonlinear=true,
                        distribution=:regular)
    DistMeshQuality(iso,
                    deltat,
                    filter_less_than,
                    termination_quality,
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
    min_volume_edge_ratio::Vector{T}
    max_volume_edge_ratio::Vector{T}
    average_volume_edge_ratio::Vector{T}
    retriangulations::Vector{Int} # Iteration num where retriangulation occured
end

DistMeshStatistics() = DistMeshStatistics{Float64}([],[],[],[],[],[])

"""
Uniform edge length function.
"""
struct HUniform end


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
