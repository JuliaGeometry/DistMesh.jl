using DistMesh
using BenchmarkTools
using Plots
using Dates
using Descartes
using GeometryBasics
#using JLD
#using HDF5

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

#
# Distance Functions
#

torus(v) = (sqrt(v[1]^2+v[2]^2)-0.5)^2 + v[3]^2 - 0.25 # torus


# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()

# all solids should be in -2:2 for each axis for ease

println("Benchmarking DistMesh.jl...")

#
# Algorithms to benchmark
#

timestamp = Dates.format(now(), "yyyy-mm-ddTHH_MM")
# generate an orthogonal test suite
for ttol in 0.01:0.01:0.05, deltat in 0.05:0.05:0.1, el = 0.1:0.05:0.2
    rt = time()
    result = distmesh(torus,
                     HUniform(),
                     el,
                     DistMeshSetup(deltat=deltat,ttol=ttol,distribution=:packed),
                     origin = GeometryBasics.Point{3,Float64}(-2),
                     widths = GeometryBasics.Point{3,Float64}(4),
                     stats=true)
    running_time = time() - rt # approximate, since we mostly care about convergence factors
    item = "torus$timestamp"
    folder = joinpath(@__DIR__, "output")
    !isdir(folder) && mkdir(folder)
    folder = joinpath(folder, "$item")
    !isdir(folder) && mkdir(folder)
    param_str = "_ttol=$(ttol)_deltat=$(deltat)_el=$(el)"
    # save plots
    plotout(result.stats, DistMesh.triangle_qualities(result.points,result.tetrahedra), folder, param_str)
    # save dataset as JLD
    # jldopen("$folder/$item.jld", "w") do file
    #     g = g_create(file, param_str)
    #     g["points"] = p
    #     g["tets"] = t
    #     g["stats"] = s
    #     g["running_time"] = running_time
    # end
    println(param_str)
end
