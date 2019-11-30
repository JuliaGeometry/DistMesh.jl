using DistMesh
using Test
using GeometryBasics

include("vals.jl")



@testset "point distributions" begin
    vlen(a,b) = sqrt(sum((a-b).^2))
    @testset "simple cubic" begin
        pts = []
        dists = []
        f(x) = -1
        DistMesh.simplecubic!(f, pts, dists, 0.5, 0, Point{3,Float64}(0),Point{3,Float64}(1),Point{3,Float64})
        @test length(pts) == 27
        @test length(dists) == 27
        @test isapprox(vlen(pts[1],pts[2]),0.5)
        @test length(pts) == length(unique(pts))
    end

    @testset "face centered cubic" begin
        pts = []
        dists = []
        f(x) = -1
        DistMesh.facecenteredcubic!(f, pts, dists, 0.5, 0, Point{3,Float64}(0),Point{3,Float64}(1),Point{3,Float64})
        @test length(pts) == 216
        @test length(dists) == 216
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

@testset "hilbert sort" begin
    rng = 0:0.3:1
    a = Vector{Vector{Float64}}(undef,length(rng)^3)
    i = 1
    for xi in rng, yi in rng, zi in rng
        a[i] = [xi,yi,zi]
        i += 1
    end
    DistMesh.hilbertsort!(a)
    @test a == hilbert_a
end

@testset "distmesh 3D" begin
    d(p) = sqrt(sum(p.^2))-1
    p,t,_ = distmesh(d,HUniform(),0.2)
    @test length(p) == 485
    @test length(t) == 2207

    p,t,_ = distmesh(d,HUniform(),0.2, DistMeshSetup(distribution=:packed))
    @test length(p) == 742
    @test length(t) == 3472

    # test stats is not messing
    p,t,s = distmesh(d,HUniform(),0.2, stats=true)
    @test length(p) == 485
    @test length(t) == 2207

    p,t,s = distmesh(d,HUniform(),0.4, stats=true)
    @test length(p) == 56
    @test length(t) == 186
    for fn in fieldnames(typeof(s))
        @test isapprox(getproperty(s,fn), getproperty(stat_04,fn))
    end
end
