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

