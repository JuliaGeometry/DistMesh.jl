
@testset "distmesh 3D" begin
    d(p) = sqrt(sum(p.^2))-1
    p,t = distmesh(d,huniform,0.2, vis=_VIS)
    @test length(p) == 485
    @test length(t) == 2207
    p,t = distmesh(d,huniform,0.2, vis=_VIS, distribution=:packed)
    @test length(p) == 742
    @test length(t) == 3455
end
