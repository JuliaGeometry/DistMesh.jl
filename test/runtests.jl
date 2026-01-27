using DistMesh
using Test
using CairoMakie

# ------------------------------------------------------------------------------
# DistMesh 2D Tests
# ------------------------------------------------------------------------------

"""
    check_mesh(p, t; np, nt, area, areatol, minqual)

Validates a 2D mesh against expected properties.
- `np`: Expected number of nodes.
- `nt`: Expected number of triangles.
- `area`: Expected total area of the domain.
- `areatol`: Tolerance for total area check.
- `minqual`: Minimum allowable triangle quality (0.0 to 1.0).
"""
function check_mesh(msh::DMesh; np::Int=0, nt::Int=0, area::Real=0.0, areatol::Real=1e-6, minqual::Real=0.1)
    @testset "Mesh Validation" begin
        
        # 1. Check Counts
        np > 0 && @test length(msh.p) == np
        nt > 0 && @test length(msh.t) == nt

        # 2. Compute Geometric Properties
        total_area = sum(element_volumes(msh))
        min_q = minimum(element_qualities(msh))
        has_negative_area = minimum(element_volumes(msh)) < 0.0

        # 3. Check Winding / Topology
        @test !has_negative_area

        # 4. Check Total Area
        area > 0.0 && @test isapprox(total_area, area, atol=areatol)

        # 5. Check Element Quality
        @test min_q >= minqual
    end
end

@testset "Unit Circle Mesh" begin
    msh = distmesh2d(dcircle, huniform, 0.2, ((-1,-1), (1,1)))
    check_mesh(msh, np=88, nt=143)
    
    for h in (0.01, 0.02, 0.05, 0.10, 0.15, 0.2)
        msh = distmesh2d(dcircle, huniform, h, ((-1,-1), (1,1)))
        check_mesh(msh,
                   area=pi, 
                   areatol=h^2, # Not an exact circle
                   minqual=0.6
                   )
    end
end


# ------------------------------------------------------------------------------
# Run all scripts in the examples/ directory

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")

@testset "Examples" begin
    # Get all .jl files
    files = filter(f -> endswith(f, ".jl"), readdir(EXAMPLES_DIR))
    
    for file in files
        @testset "$file" begin
            # Read the example file
            path = joinpath(EXAMPLES_DIR, file)
            content = read(path, String)
            
            # Swap GLMakie -> CairoMakie (for Headless tests)
            content = replace(content, "using GLMakie" => "using CairoMakie")

            # Run the code
            try
                include_string(Main, content, path)
                @test true 
            catch e
                @error "Example $file failed to run" exception=(e, catch_backtrace())
                @test false
            end
        end
    end
end
