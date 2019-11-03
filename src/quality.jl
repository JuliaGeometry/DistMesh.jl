"""
Determine the quality of a triangle given 3 points.

Points must be 3D.
"""
function triqual(p1, p2, p3)
    d12 = p2 - p1
    d13 = p3 - p1
    d23 = p3 - p2
    n   = cross(d12,d13)
    vol = sqrt(sum(n.^2))/2
    den = dot(d12,d12) + dot(d13,d13) + dot(d23,d23)
    return sqrt(3)*4*vol/den
end

function triangle_qualities(p,tets)
    tris = NTuple{3,Int}[]
    tets_to_tris!(tris, tets)
    qualities = Vector{Float64}(undef,length(tris))
    triangle_qualities!(tris,qualities,p,tets)
end

function triangle_qualities!(tris,triset,qualities,p,tets)
    tets_to_tris!(tris,triset,tets)
    resize!(qualities, length(tris))
    for i in eachindex(tris)
        tp = tris[i]
        qualities[i] = triqual(p[tp[1]], p[tp[2]], p[tp[3]])
    end
    qualities
end

function triangle_qualities!(tris::Vector,qualities::Vector,p,tets)
    triangle_qualities!(tris,Set{eltype(tris)}(),qualities,p,tets)
end