module DistMesh

using Combinatorics,
      GeometryTypes,
      Meshing,
      LinearAlgebra,
      StaticArrays,
      TetGen

#include("circumcenter.jl")
#include("distmeshsurface.jl")
include("huniform.jl")
include("dsphere.jl")
include("mkt2t.jl")
include("munique.jl")
include("distmeshnd.jl")
#include("trisurfupd.jl")

export distmeshsurface
export distmeshnd

end # module
