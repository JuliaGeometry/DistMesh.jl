
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
                push!(triset, bitpack(p1,p2,p3))
            elseif p1 <= p3 <= p2
                push!(triset, bitpack(p1,p3,p2))
            elseif p2 <= p3 <= p1
                push!(triset, bitpack(p2,p3,p1))
            elseif p2 <= p1 <= p3
                push!(triset, bitpack(p2,p1,p3))
            elseif p3 <= p1 <= p2
                push!(triset, bitpack(p3,p1,p2))
            elseif p3 <= p2 <= p1
                push!(triset, bitpack(p3,p2,p1))
            end
        end
    end
    resize!(tris, length(triset))
    #copy elements to tri
    i = 1
    for elt in triset
        tris[i] = unbitpack3(elt)
        i = i + 1
    end
    sort!(tris)
    tris
end

function tets_to_tris!(tris::Vector, tets::Vector)
    tets_to_tris!(tris, Set{UInt64}(), tets)
end

"""
    Decompose tets to edges, using a pre-allocated array and set.
    Set ensures uniqueness, and result will be sorted.
"""
function tet_to_edges!(pair::Vector, pair_set::Set, t)
    empty!(pair_set)
    for i in eachindex(t)
        for ep in 1:6
            p1 = t[i][tetpairs[ep][1]]
            p2 = t[i][tetpairs[ep][2]]
            push!(pair_set, p1 > p2 ? bitpack(p2,p1) : bitpack(p1,p2))
        end
    end
    resize!(pair, length(pair_set))
    # copy pair set to array since sets are not sortable
    i = 1
    for elt in pair_set
        pair[i] = unbitpack2(elt)
        i = i + 1
    end

    # sort the edge pairs for better point lookup
    sort!(pair)
end

"""
Pack two 32 bit numbers (or smaller) in an UInt64 to speed up hashing in dicts
"""
function bitpack(xi,yi)
    (unsafe_trunc(UInt64, xi) << 32) | unsafe_trunc(UInt64, yi)
end

function bitpack(xi,yi,zi)
    (unsafe_trunc(UInt64, xi) << 42) | (unsafe_trunc(UInt64, yi & 0x01fffff) << 21) | unsafe_trunc(UInt64, zi & 0x01fffff)
end

function unbitpack2(key)
    # we can use unsafe truncation here because the bounds
    # are enforced in the initial masking.
    unsafe_trunc(UInt32, key >>> 32), unsafe_trunc(UInt32, key)
end

function unbitpack3(key)
    # we can use unsafe truncation here because the bounds
    # are enforced in the initial masking.
    unsafe_trunc(UInt32, key >>> 42), unsafe_trunc(UInt32, key >>> 21) & 0x01fffff, unsafe_trunc(UInt32, key) & 0x01fffff
end