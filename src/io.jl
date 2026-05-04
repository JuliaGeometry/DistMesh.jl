"""
    read_stl(fname)

Read mesh from STL (Stereolithography) file.

Reads vertices and faces from an STL file and returns a DMesh.
Currently only supports ASCII STL format.

# Arguments
- `fname`: Input filename

# Returns
- `DMesh`: Mesh with vertices and elements

# Example
```julia
msh = read_stl("model.stl")
```
"""
function read_stl(fname)
    p = Point3d[]
    for line in eachline(fname)
        line = lowercase(strip(line))
        if startswith(line, "vertex")
            parts = split(line)
            length(parts) == 4 || error("Syntax error: Malformed vertex definition.")
            push!(p, Tuple(parse.(Float64, @view parts[2:4])))
        end
    end
    t = [ Index3(i, i+1, i+2) for i in 1:3:length(p) ]
    return DMesh(p, t)
end

"""
    write_ply(fname, msh::DMesh{D,T,E,N}) where {D,T,E,N}

Export mesh to PLY (Polygon File Format) file.

The PLY format is written in ASCII 1.0 format with vertices and faces. 
Vertex indices are converted from Julia's 1-based indexing to PLY's 0-based indexing.

# Arguments
- `msh::DMesh{D,T,E,N}`: Mesh to export
- `fname`: Output filename

# Example
```julia
write_ply("output.ply", msh)
```
"""
function write_ply(fname, msh::DMesh{D,T,E,N}) where {D,T,E,N}
    nv = length(msh.p)
    nf = length(msh.t)
    
    # Dimension names for property headers
    dim_names = ["x", "y", "z", "w"]
    
    open(fname, "w") do io
        # Write header
        println(io, "ply")
        println(io, "format ascii 1.0")
        println(io, "element vertex ", nv)
        
        # Write vertex properties based on dimension
        for i in 1:D
            println(io, "property float ", dim_names[i])
        end
        
        println(io, "element face ", nf)
        println(io, "property list uchar int vertex_index")
        println(io, "end_header")
        
        # Write vertices
        for p in msh.p
            println(io, join(p, " "))
        end
        
        # Write faces (convert to 0-based indexing)
        for t in msh.t
            println(io, N, " ", join(t .- 1, " "))
        end
    end
end

"""
    read_ply(fname)

Read mesh from PLY (Polygon File Format) file.

Reads vertices and faces from a PLY file and returns a DMesh.
Vertex indices are converted from PLY's 0-based indexing to Julia's 1-based indexing.
The mesh dimension and element type are inferred from the file header and first face.

# Arguments
- `fname`: Input filename

# Returns
- `DMesh`: Mesh with vertices and elements

# Example
```julia
msh = read_ply("input.ply")
```
"""
function read_ply(fname)
    open(fname, "r") do io
        # Parse header to determine dimension and element count
        dim = 0
        nv = 0
        
        for line in eachline(io)
            line = strip(line)
            line == "end_header" && break
            
            if startswith(line, "element vertex")
                nv = parse(Int, split(line)[3])
            elseif occursin(r"^property (float|double) ([xyzw])$", line)
                dim += 1
            end
        end
        
        # Read vertices with proper typing
        p_data = Vector{NTuple{dim, Float64}}(undef, nv)
        for i in 1:nv
            line = readline(io)
            coords = parse.(Float64, split(strip(line)))
            p_data[i] = tuple(coords[1:dim]...)
        end
        
        # Read first face to determine element size N
        first_face_line = readline(io)
        parts = parse.(Int, split(strip(first_face_line)))
        n = parts[1]
        
        # Create properly typed vectors
        p = Vector{NTuple{dim, Float64}}(undef, nv)
        p .= p_data
        
        t = Vector{NTuple{n, Int}}(undef, 0)
        sizehint!(t, nv)  # Approximate guess on face count
        
        # Add first face
        indices = (@view parts[2:n+1]) .+ 1  # Convert 0-based to 1-based
        push!(t, tuple(indices...))
        
        # Read remaining faces
        for line in eachline(io)
            line = strip(line)
            isempty(line) && continue
            
            parts = parse.(Int, split(line))
            n = parts[1]
            indices = (@view parts[2:n+1]) .+ 1  # Convert 0-based to 1-based
            push!(t, tuple(indices...))
        end
        
        return DMesh(p, t)
    end
end

"""
    write_stl(fname, msh::DMesh{3,T,Simplex{2},3}; meshname="DistMesh_Exported_Mesh")

Export 3D triangular mesh to STL (Stereolithography) ASCII format.

Writes vertices and triangular faces to an STL file. Surface normals are automatically 
calculated from the three vertices of each triangle using the cross product.

# Arguments
- `msh::DMesh{3,T,Simplex{2},3}`: 3D triangular mesh to export
- `fname`: Output filename
- `meshname::String`: Name for the solid object (default: "DistMesh_Exported_Mesh")

# Example
```julia
write_stl("output.stl", msh, meshname="MyModel")
```
"""
function write_stl(fname, msh::DMesh{3,T,Simplex{2},3}; meshname="DistMesh_Exported_Mesh") where T
    open(fname, "w") do io
        println(io, "solid \"$meshname\"")
        
        for tri in msh.t
            # Get the three vertices
            p1 = msh.p[tri[1]]
            p2 = msh.p[tri[2]]
            p3 = msh.p[tri[3]]
            
            # Calculate normal from cross product
            v1 = p2 .- p1
            v2 = p3 .- p1
            normal = cross(v1, v2)
            normal_normalized = normal ./ norm(normal)
            
            # Write facet
            println(io, "  facet normal ", join(normal_normalized, " "))
            println(io, "    outer loop")
            println(io, "      vertex ", join(p1, " "))
            println(io, "      vertex ", join(p2, " "))
            println(io, "      vertex ", join(p3, " "))
            println(io, "    endloop")
            println(io, "  endfacet")
        end
        
        println(io, "endsolid \"$meshname\"")
    end
end

