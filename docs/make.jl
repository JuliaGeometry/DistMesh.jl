using Documenter, DistMesh

makedocs(
    modules=[DistMesh],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    sitename="DistMesh.jl",
    authors="Steve Kelly <kd2cca@gmail.com>",
    assets=String[],
)

deploydocs(
    repo="github.com/JuliaGeometry/DistMesh.jl",
)
