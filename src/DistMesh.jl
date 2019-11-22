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
    ttol::T
    ptol::T
    distribution::Symbol # intial point distribution
end

function DistMeshSetup(;iso=0,
                        ptol=.001,
                        deltat=0.05,
                        ttol=0.02,
                        distribution=:regular)
    T = promote_type(typeof(iso),typeof(ptol),typeof(deltat), typeof(ttol))
    DistMeshSetup{T}(iso,
                     deltat,
                     ttol,
                     ptol,
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
    retriangulations::Vector{Int} # Iteration num where retriangulation occured
end

DistMeshStatistics() = DistMeshStatistics{Float64}([],[],[],[],[])

"""
Uniform edge length function.
"""
struct HUniform end

#include("circumcenter.jl")
#include("distmeshsurface.jl")
include("diff.jl")
include("pointdistribution.jl")
include("mkt2t.jl")
include("distmeshnd.jl")
include("compat/munique.jl")
include("tetgen.jl")
include("quality.jl")
include("decompositions.jl")
#include("trisurfupd.jl")

#export distmeshsurface
export distmesh, DistMeshSetup, DistMeshStatistics, HUniform

end # module
