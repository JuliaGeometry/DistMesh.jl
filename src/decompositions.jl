
"""
convert tets to tris,
returned sorted and unique
"""
function tets_to_tris!(tris, triset, tets)
    empty!(triset)
    for i in eachindex(tets)
        for j in 1:4
            tp = tettriangles[j]
            p1 = tets[i][tp[1]]
            p2 = tets[i][tp[2]]
            p3 = tets[i][tp[3]]
            # sort indices
            if p1 <= p2 <= p3
                push!(triset, (p1,p2,p3))
            elseif p1 <= p3 <= p2
                push!(triset, (p1,p3,p2))
            elseif p2 <= p3 <= p1
                push!(triset, (p2,p3,p1))
            elseif p2 <= p1 <= p3
                push!(triset, (p2,p1,p3))
            elseif p3 <= p1 <= p2
                push!(triset, (p3,p1,p2))
            elseif p3 <= p2 <= p1
                push!(triset, (p3,p2,p1))
            end
        end
    end
    resize!(tris, length(triset))
    #copy elements to tri
    i = 1
    for elt in triset
        tris[i] = elt
        i = i + 1
    end
    sort!(tris)
    tris
end

function tets_to_tris!(tris::Vector, tets::Vector)
    tets_to_tris!(tris, Set{eltype(tris)}(), tets)
end