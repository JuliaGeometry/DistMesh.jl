module DistMesh

using Meshing,
      LinearAlgebra,
      StaticArrays,
      TetGen,
      Makie,
      AbstractPlotting,
      Colors
using GeometryBasics
using GeometryBasics: Triangle, Tetrahedron, Mesh, Polytope, Point
import GeometryTypes

_DEFAULT_SAMPLES = (24,24,24)
const tetpairs = ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))

struct DistMeshSetup{T}
    iso::T
    deltat::T
    ttol::T
    ptol::T
end

#include("circumcenter.jl")
#include("distmeshsurface.jl")
include("huniform.jl")
include("mkt2t.jl")
include("distmeshnd.jl")
include("compat/munique.jl")
include("compat/delaunayn.jl")
#include("trisurfupd.jl")

#export distmeshsurface
export distmeshnd

export huniform

end # module
