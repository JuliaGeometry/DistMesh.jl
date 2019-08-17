module DistMesh

using GeometryTypes,
      Meshing,
      StaticArrays
include("GroupSlices.jl")
using .GroupSlices

#include("circumcenter.jl")
#include("distmeshsurface.jl")
include("dsphere.jl")
include("mkt2t.jl")

end # module
