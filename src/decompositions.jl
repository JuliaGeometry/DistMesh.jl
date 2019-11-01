
"""
convert tets to tris,
returned sorted and unique
"""
function tets_to_tris!(tris, tets)
    resize!(tris, length(tets)*4)
    for i in eachindex(tets)
        for j in 1:4
            tp = tettriangles[j]
            p1 = tets[i][tp[1]]
            p2 = tets[i][tp[2]]
            p3 = tets[i][tp[3]]
            # sort indices
            if p1 <= p2 <= p3
                tris[(i-1)*4+j] = (p1,p2,p3)
            elseif p1 <= p3 <= p2
                tris[(i-1)*4+j] = (p1,p3,p2)
            elseif p2 <= p3 <= p1
                tris[(i-1)*4+j] = (p2,p3,p1)
            elseif p2 <= p1 <= p3
                tris[(i-1)*4+j] = (p2,p1,p3)
            elseif p3 <= p1 <= p2
                tris[(i-1)*4+j] = (p3,p1,p2)
            elseif p3 <= p2 <= p1
                tris[(i-1)*4+j] = (p3,p2,p1)
            end
        end
    end
    sort!(tris)
    unique!(tris)
    tris
end