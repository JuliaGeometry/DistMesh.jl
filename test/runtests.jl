using DistMesh
using Test
using GeometryBasics

include("vals.jl")

const DistMeshND = DistMesh.DistMeshND

@testset "point distributions" begin
    vlen(a,b) = sqrt(sum((a-b).^2))
    @testset "simple cubic" begin
        pts = []
        dists = []
        f(x) = -1
        DistMeshND.simplecubic!(f, pts, dists, 0.5, 0, Point{3,Float64}(0),Point{3,Float64}(1),Point{3,Float64})
        @test length(pts) == 27
        @test length(dists) == 27
        @test isapprox(vlen(pts[1],pts[2]),0.5)
        @test length(pts) == length(unique(pts))
    end

    @testset "face centered cubic" begin
        pts = []
        dists = []
        f(x) = -1
        DistMeshND.facecenteredcubic!(f, pts, dists, 0.5, 0, Point{3,Float64}(0),Point{3,Float64}(1),Point{3,Float64})
        @test length(pts) == 216
        @test length(dists) == 216
        @test isapprox(vlen(pts[1],pts[2]),0.5)
        @test length(pts) == length(unique(pts))
    end

end

@testset "quality analysis" begin
    @testset "triangles" begin
        @test DistMeshND.triqual([0,0,0],[1,0,0],[0,1,0]) == DistMeshND.triqual([0,0,0],[2,0,0],[0,2,0])
        @test DistMeshND.triqual([0,0,0],[1,0,1],[0,1,1]) == DistMeshND.triqual([0,0,0],[2,0,2],[0,2,2])
        @test DistMeshND.triqual([0,0,0],[2,0,0],[1,sqrt(3),0]) ≈ 1
        @test DistMeshND.triqual([0,0,0],[1,sqrt(3),0],[2,0,0]) ≈ 1
    end
    @testset "volume-length" begin
        pts = ([1,0,-1/sqrt(2)], [-1,0,-1/sqrt(2)], [0,1,1/sqrt(2)], [0,-1,1/sqrt(2)])
        pts2 = ([1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1])
        pts_degenerate = ([1,1,1], [1,1,1], [-1,1,-1], [-1,-1,1])
        @test DistMeshND.volume_edge_ratio(pts...) ≈ 1
        @test DistMeshND.volume_edge_ratio((pts.*2)...) ≈ 1
        @test DistMeshND.volume_edge_ratio((pts.*1e-6)...) ≈ 1
        @test isnan(DistMeshND.volume_edge_ratio((pts.*0)...))
        @test DistMeshND.volume_edge_ratio(pts2...) ≈ 1
        @test DistMeshND.volume_edge_ratio((pts2.*2)...) ≈ 1
        @test DistMeshND.volume_edge_ratio(pts_degenerate...) == 0

    end
end

@testset "decompositions" begin
    @testset "tets to triangles" begin
        simps = [[4,3,2,1],[5,4,3,2],[1,2,3,4]]
        tris = Tuple{Int,Int,Int}[]
        DistMeshND.tets_to_tris!(tris,simps)
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
    DistMeshND.hilbertsort!(a)
    @test a == hilbert_a
end

@testset "distmesh 3D" begin
    d(p) = sqrt(sum(p.^2))-1
    result = distmeshnd(d,HUniform(),0.2)
    @test length(result.points) == 485
    @test length(result.tetrahedra) == 2207

    result = distmeshnd(d,HUniform(),0.2, DistMeshSetup(distribution=:packed))
    @test length(result.points) == 742
    @test length(result.tetrahedra) == 3472

    # test stats is not messing
    result = distmeshnd(d,HUniform(),0.2, stats=true)
    @test length(result.points) == 485
    @test length(result.tetrahedra) == 2207

    result = distmeshnd(d,HUniform(),0.4, stats=true)
    @test length(result.points) == 56
    @test length(result.tetrahedra) == 186
    #for fn in fieldnames(typeof(result.stats))
    #    @test isapprox(getproperty(result.stats,fn), getproperty(stat_04,fn))
    #end
end

@testset "dihedral metrics" begin
    d(p) = sqrt(sum(p.^2))-1
    result = distmeshnd(d,HUniform(),0.2)
    p = result.points
    t = result.tetrahedra
    all_angs =  DistMeshND.dihedral_angles(p,t)
    min_angs =  DistMeshND.min_dihedral_angles(p,t)
    ax = extrema(all_angs)
    mx = extrema(min_angs)
    @test ax[1] == mx[1]
    @test all(ax .≈ (0.023502688273828173, 3.104396619996953))
    @test all(mx .≈ (0.023502688273828173, 1.2044902180168893))
end
