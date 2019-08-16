using Documenter, DistMesh

makedocs(;
    modules=[DistMesh],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/sjkelly/DistMesh.jl/blob/{commit}{path}#L{line}",
    sitename="DistMesh.jl",
    authors="Steve Kelly <kd2cca@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/sjkelly/DistMesh.jl",
)
