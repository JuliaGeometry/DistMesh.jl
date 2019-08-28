module DistMesh

using GeometryTypes,
      Meshing,
      StaticArrays

#include("circumcenter.jl")
#include("distmeshsurface.jl")
include("dsphere.jl")
include("mkt2t.jl")
include("munique.jl")
include("trisurfupd.jl")

export distmeshsurface

end # module
