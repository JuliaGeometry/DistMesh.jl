using DistMesh
using BenchmarkTools
using Plots
using Dates
using Descartes

include("util.jl")

#=
In addition to performance we want to measure the convergence rate.
Ideally this will be somewhat smooth without many fast changes.

- ttol (retriangulation critera)
- deltat (displacement applied per iteration)

termination critera:
- ptol (termination critera)
- mamove delta (TBD)

output for each trial:
- store input parameters (above)
- DistMeshStatistics
- Benchmark Stats
=#

#=
TODO:

- Make sure plots have equal yscale
- Consistent quality histogram bins

Solids:
- Box
- Sphere
- Shelled Sphere
- Torus
- deadmau5 (mouse-head with sphere subtracted from sphere, bad metric distance function)
=#

# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()

# all solids should be in -2:2 for each axis for ease
suite["Sphere"] = BenchmarkGroup()
suite["Torus"] = BenchmarkGroup()

println("Benchmarking DistMesh.jl...")

#
# Algorithms to benchmark
#

algos = []

# generate an orthogonal test suite
for ttol in 0.1:0.2:0.2

end

#
# Benchmark SDF constructon for SDF and Direct sample comparisons
#

suite["SDF Construction"]["Torus"] = @benchmarkable SignedDistanceField(HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)),0.05) do v
    (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus
end

#
# Solid Constructions
#

torus(v) = (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus

#
# Benchmark algorithms
#

for algo in algos_sdf
    suite["SDF Mesh"][string(typeof(algo))] = @benchmarkable HomogenousMesh(sdf_torus, $algo)
end

for algo in algos_fn
    suite["Function Mesh"][string(typeof(algo))] = @benchmarkable HomogenousMesh(fn_torus, HyperRectangle(Vec(-2,-2,-2.),Vec(4,4,4.)), (81,81,81), $algo)
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

