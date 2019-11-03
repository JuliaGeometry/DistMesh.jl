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
struct DistMeshSetup{T,RT}
    iso::T
    deltat::T
    retriangulation_criteria::RT
    ptol::T
    droptets::Bool # drop tetrahedra with centroids outside of the boundary
end

function DistMeshSetup(;iso=0,
                        ptol=.001,
                        deltat=0.05,
                        retriangulation_criteria=RetriangulateMaxMove(0.02),
                        droptets=true)
    T = promote_type(typeof(iso),typeof(ptol),typeof(deltat))
    DistMeshSetup{T,typeof(retriangulation_criteria)}(iso,
                                                      deltat,
                                                      retriangulation_criteria,
                                                      ptol,
                                                      droptets)
end

abstract type AbstractRetriangulationCriteria end

"""
    retriangulate if a move over a given tolerance occurs
"""
struct RetriangulateMaxMove{T} <: AbstractRetriangulationCriteria
    ttol::T
end

"""
    retriangulate on a positive delta over N moves, after M force iterations
"""
struct RetriangulateMaxMoveDelta <: AbstractRetriangulationCriteria
    move_count::Int
    iterations::Int
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

export huniform

end # module
