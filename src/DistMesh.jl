module DistMesh

using Meshing,
      LinearAlgebra,
      StaticArrays,
      TetGen
using GeometryBasics
using GeometryBasics: Triangle, Tetrahedron, Mesh, Polytope, Point
import GeometryTypes

_DEFAULT_SAMPLES = (24,24,24)
const tetpairs = ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
const tettriangles = ((1,2,3),(1,2,4),(2,3,4),(1,3,4))


"""
    DistMeshSetup

    iso (default: 0): Value of which to extract the iso surface, inside negative
    deltat (default: 0.1): the fraction of edge displacement to apply each iteration
"""
struct DistMeshSetup{T}
    iso::T
    deltat::T
    ptol::T
    ttol::T
    maxmove_delta::Int
    maxmove_delta_delay::Int
    droptets::Bool # drop tetrahedra with centroids outside of the boundary
    initial_points::Symbol
end

function DistMeshSetup(;iso=0.0,
                        ptol=.001,
                        deltat=0.2,
                        ttol=0.1,
                        maxmove_delta=6, # number of prior samples to consider
                        maxmove_delta_delay=2,
                        droptets=true,
                        initial_points=:regular)
    T = promote_type(typeof(iso),typeof(ptol),typeof(deltat))
    DistMeshSetup{T}(iso,
                     deltat,
                     ptol,
                     ttol,
                     maxmove_delta,
                     maxmove_delta_delay,
                     droptets,
                     initial_points)
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
    retriangulations::Vector{Int} # Iteration num where retriangulation occured
end

DistMeshStatistics() = DistMeshStatistics{Float64}([],[],[],[],[])


#include("circumcenter.jl")
#include("distmeshsurface.jl")
include("pointdistribution.jl")
include("huniform.jl")
include("mkt2t.jl")
include("distmeshnd.jl")
include("compat/munique.jl")
include("tetgen.jl")
include("quality.jl")
include("decompositions.jl")
#include("trisurfupd.jl")

#export distmeshsurface
export distmesh, DistMeshSetup, DistMeshStatistics
export RetriangulateMaxMove, RetriangulateMaxMoveDelta
export huniform

end # module
