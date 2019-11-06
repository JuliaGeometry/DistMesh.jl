using DistMesh
using Test
using MAT
using GeometryTypes
using StaticArrays


@testset "distmesh 3D" begin
    d(p) = sqrt(sum(p.^2))-1
    p,t,_ = distmesh(d,huniform,0.2)
    @test length(p) == 485
    @test length(t) == 2207

    p,t,_ = distmesh(d,huniform,0.2, DistMeshSetup(initial_points=:packed))
    @test length(p) == 742
    @test length(t) == 3472
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


@testset "point distributions" begin
    vlen(a,b) = sqrt(sum((a-b).^2))
    @testset "simple cubic" begin
        pts = []
        f(x) = -1
        DistMesh.simplecubic!(f, pts, 0.5, 0, Point{3,Float64}(0),Point{3,Float64}(1),Point{3,Float64})
        @test length(pts) == 27
        @test isapprox(vlen(pts[1],pts[2]),0.5)
        @test length(pts) == length(unique(pts))
    end

    @testset "face centered cubic" begin
    pts = []
    f(x) = -1
    DistMesh.facecenteredcubic!(f, pts, 0.5, 0, Point{3,Float64}(0),Point{3,Float64}(1),Point{3,Float64})
    @test length(pts) == 216
    @test isapprox(vlen(pts[1],pts[2]),0.5)
    @test length(pts) == length(unique(pts))
    end

end

@testset "quality analysis" begin
    @testset "triangles" begin
        @test DistMesh.triqual([0,0,0],[1,0,0],[0,1,0]) == DistMesh.triqual([0,0,0],[2,0,0],[0,2,0])
        @test DistMesh.triqual([0,0,0],[1,0,1],[0,1,1]) == DistMesh.triqual([0,0,0],[2,0,2],[0,2,2])
        @test DistMesh.triqual([0,0,0],[2,0,0],[1,sqrt(3),0]) ≈ 1
        @test DistMesh.triqual([0,0,0],[1,sqrt(3),0],[2,0,0]) ≈ 1
    end

end

@testset "decompositions" begin
    @testset "tets to triangles" begin
        simps = [[4,3,2,1],[5,4,3,2],[1,2,3,4]]
        tris = Tuple{Int,Int,Int}[]
        DistMesh.tets_to_tris!(tris,simps)
        @test tris == [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4), (2, 3, 5), (2, 4, 5), (3, 4, 5)]
    end
end

# @testset "distmeshsurface" begin
#     fd(p) = dsphere(p,0,0,0,1);
#     fh(p) = 0.05+0.5*dsphere(p,0,0,1,0);
#     p,t = distmeshsurface(fd, fh, 0.15, 1.1.*[-1 -1 -1;1 1 1]);
#     @show p, t
# end
