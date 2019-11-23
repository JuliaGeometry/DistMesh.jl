
function simplecubic!(fdist, points, dists, h, iso, origin, widths, ::Type{VertType}) where VertType
    @inbounds for xi = origin[1]:h:(origin[1]+widths[1]), yi = origin[2]:h:(origin[2]+widths[2]), zi = origin[3]:h:(origin[3]+widths[3])
        point = VertType(xi,yi,zi)
        d = fdist(point)
        if d < iso
            push!(points,point)
            push!(dists, d)
        end
    end
end

function facecenteredcubic!(fdist, points, dists, h, iso, origin, widths, ::Type{VertType}) where VertType
    # face-centered cubic point distribution
    r = h/2
    counts = round.(widths./h).+2
    @inbounds for xi = -1:Int(counts[1]), yi = -1:Int(counts[2]), zi = -1:Int(counts[3])
        point = VertType(2xi+((yi+zi)%2), sqrt(3)*(yi+(zi%2)/3),2*sqrt(6)*zi/3).*r + origin
        d = fdist(point)
        if d < iso
            push!(points,point)
            push!(dists, d)
        end
    end
end