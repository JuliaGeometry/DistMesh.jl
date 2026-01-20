using Documenter
using DistMesh
using Literate

# --- 1. Automation Script: Build the "Examples" Page ---

# Define paths
examples_dir = joinpath(@__DIR__, "..", "examples")
docs_src     = joinpath(@__DIR__, "src")
master_file  = joinpath(docs_src, "examples.md")

# Initialize the master "Examples" page
open(master_file, "w") do io
    write(io, "# Examples\n\nA collection of 2D meshing examples.\n\n")
end

# Find all .jl files and sort them
jl_files = sort(filter(f -> endswith(f, ".jl"), readdir(examples_dir)))

for file in jl_files
    source_path = joinpath(examples_dir, file)
    dest_path   = joinpath(docs_src, file) # Copy directly to docs/src/
    
    # 1. Copy the .jl file to docs/src/
    #    This satisfies Documenter, which looks for the source file next to examples.md
    cp(source_path, dest_path, force=true)
    
    # 2. Run Literate
    #    We output the .md file to the same folder (docs/src/) to keep paths simple.
    #    We use name=... to ensure unique scope.
    Literate.markdown(dest_path, docs_src; 
        documenter=true, 
        credit=false,
        name=replace(file, ".jl"=>"")
    )
    
    # 3. Read the generated markdown
    md_file = replace(file, ".jl" => ".md")
    generated_md_path = joinpath(docs_src, md_file)
    content = read(generated_md_path, String)
    
    # 4. Append it to the master examples.md file
    open(master_file, "a") do io
        write(io, "\n\n" * content)
    end

    # 5. Cleanup the individual .md file (we don't need it anymore, it's in examples.md)
    #    Note: We KEEP the .jl file in docs/src/ so the "Source" links work.
    rm(generated_md_path)
end

# --- 2. Build Documentation ---

makedocs(
    sitename = "DistMesh.jl",
    modules  = [DistMesh],
    authors  = "Steve Kelly <kd2cca@gmail.com> and Per-Olof Persson <persson@berkeley.edu>",
    pages    = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
        "Legacy (N-D & Theory)" => "distmeshnd.md",
    ],
    warnonly = [:missing_docs],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/JuliaGeometry/DistMesh.jl.git",
    push_preview = true
)

