
@testset begin
    d(p) = sqrt(sum(p.^2))-1
    p,t = distmeshnd(d,huniform,0.2)
    @test length(p) == 485
    @test length(t) == 2057
end
