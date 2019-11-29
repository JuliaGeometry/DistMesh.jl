using DistMesh
using BenchmarkTools
using StaticArrays

a = [rand(SVector{3,Float64}) for i = 1:50000]

function test_hilbert(a)
    b = copy(a)
    DistMesh.hilbertsort!(b)
end

@benchmark test_hilbert(a)
