


@testset "mkt2t" begin
    # load up matlab workspace for comparison
    # based on sphere case
    vars = matread("mats/mkt2t.mat")
    test_t2t, test_t2n = DistMesh.mkt2t(vars["t"])
    @test test_t2n == vars["t2n"]
    @test test_t2t == vars["t2t"]

    # close?
end