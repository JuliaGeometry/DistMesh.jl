module DistMesh

using Combinatorics,
      Meshing,
      LinearAlgebra,
      StaticArrays,
      TetGen
using GeometryBasics: Triangle, Tetrahedron, Mesh, Polytope, Point
import GeometryTypes

_DEFAULT_SAMPLES = (24,24,24)

#include("circumcenter.jl")
#include("distmeshsurface.jl")
include("huniform.jl")
include("dsphere.jl")
include("mkt2t.jl")
include("distmeshnd.jl")
include("compat/munique.jl")
include("compat/delaunayn.jl")
#include("trisurfupd.jl")

#export distmeshsurface
export distmeshnd

export huniform

end # module
