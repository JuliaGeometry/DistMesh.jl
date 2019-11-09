using DistMesh
using BenchmarkTools
using GeometryBasics

#
# DistMesh.jl Performance Benchmarks
#

# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()
suite["Torus"] = BenchmarkGroup()
suite["Sphere"] = BenchmarkGroup()
suite["Box"] = BenchmarkGroup()

println("Benchmarking DistMesh.jl...")

#
# Algorithms to benchmark
#

algos = [DistMeshSetup(deltat=0.1, distribution=:packed),DistMeshSetup(distribution=:packed)]

fn_torus(v) = (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus
fn_sphere(v) = sqrt(sum(v.^2)) -1
#
# Benchmark algorithms
#

for algo in algos
    for el in (0.15,0.2)
        suite["Torus"][string(algo)*" edge="*string(el)] =
            @benchmarkable distmesh(fn_torus,
                                    huniform,
                                    $el,
                                    $algo,
                                    origin = GeometryBasics.Point{3,Float64}(-2),
                                    widths = GeometryBasics.Point{3,Float64}(4),
                                    stats=false)
        suite["Sphere"][string(algo)*" edge="*string(el)] =
            @benchmarkable distmesh(fn_sphere,
                                    huniform,
                                    $el,
                                    $algo,
                                    origin = GeometryBasics.Point{3,Float64}(-1),
                                    widths = GeometryBasics.Point{3,Float64}(2),
                                    stats=false)
    end
end

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals);
else
    tune!(suite)
    BenchmarkTools.save(paramspath, params(suite));
end

#
# Perform benchmarks and print results
#

results = run(suite)

for trial in results
    ctx = IOContext(stdout, :verbose => true, :compact => false)
    println(ctx)
    println(ctx, trial.first)
    println(ctx, trial.second)
end
