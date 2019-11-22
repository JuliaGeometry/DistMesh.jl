using DistMesh
using GeometryBasics
using PProf

fn_sphere(v) = sqrt(sum(v.^2)) -1
algo = DistMeshSetup(distribution=:packed)

distmesh(fn_sphere,HUniform(),0.15,algo,origin = GeometryBasics.Point{3,Float64}(-1),widths = GeometryBasics.Point{3,Float64}(2),stats=false)

@pprof distmesh(fn_sphere,HUniform(),0.15,algo,origin = GeometryBasics.Point{3,Float64}(-1),widths = GeometryBasics.Point{3,Float64}(2),stats=false)
