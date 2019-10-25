
@testset begin
    d(p) = sqrt(sum(p.^2))-1
    p,t = distmeshnd(d,huniform,0.2, vis=_VIS)
    @test length(p) == 485
    @test length(t) == 2207
end
