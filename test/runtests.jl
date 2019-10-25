using DistMesh
using Test
using MAT

# use makie to visualize triangulations
_VIS = false


@testset "distmesh 3D" begin
    d(p) = sqrt(sum(p.^2))-1
    p,t = distmesh(d,huniform,0.2, vis=_VIS)
    @test length(p) == 485
    @test length(t) == 2207

    p,t = distmesh(d,huniform,0.2, vis=_VIS, distribution=:packed)
    @test length(p) == 742
    @test length(t) == 3455
end


# translation utils
@testset "munique" begin
    a = [1 2; 3 4; 5 6; 1 2; 3 4; 8 8]
    @test DistMesh.munique(a) == [1,2,3,1,2,4]
end

# distmesh utils
@testset "mkt2t" begin
    # load up matlab workspace for comparison
    # based on sphere case
    vars = matread((@__DIR__)*"/mats/mkt2t.mat")
    test_t2t, test_t2n = DistMesh.mkt2t(vars["t"])
    @test test_t2n == vars["t2n"]
    @test test_t2t == vars["t2t"]

    # close?
end



# @testset "distmeshsurface" begin
#     fd(p) = dsphere(p,0,0,0,1);
#     fh(p) = 0.05+0.5*dsphere(p,0,0,1,0);
#     p,t = distmeshsurface(fd, fh, 0.15, 1.1.*[-1 -1 -1;1 1 1]);
#     @show p, t
# end
