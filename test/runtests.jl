using DistMesh
using Test
using FileIO
using MAT
using WriteVTK
using Makie

include("distmeshnd.jl")

# translation utils
include("munique.jl")

# distmesh utils
include("mkt2t.jl")
include("trisurfupd.jl")

# @testset "distmeshsurface" begin
#     fd(p) = dsphere(p,0,0,0,1);
#     fh(p) = 0.05+0.5*dsphere(p,0,0,1,0);
#     p,t = distmeshsurface(fd, fh, 0.15, 1.1.*[-1 -1 -1;1 1 1]);
#     @show p, t
# end
