
"""
A set implementation optimized for performance and integer tuples.
Uses a lot of memory, so on most computers this will handle fewer than 15000 elements for a two-tuple
"""
struct TupleSet{T} where {T}
    slots::Vector{UInt16} # stores the index of the element
    vals::Vector{T}
    bitsize::UInt8
end

"""
Intialize a tuple set where the maximum element is of `n` length.
"""
function TupleSet{T}(n::Int) where {T}
    bitsize = ceil(log2(n^length(T.parameters)))
    bitsize > 25 && error("TupleSet too large")
    TupleSet(fill(zero(UInt32), 2^bitsize), T[], bitsize)
end

function tuplehash(t::Tuple, sz)
    res = UInt16(1)
    for i in eachindex(t)
        res = sz*res + UInt16(t[i])
    end
    res
end

function push!(set::TupleSet{T}, val::T) where {T}
    ind = tuplehash(val, set.bitsize)
    if iszero(set.slots[ind]) # we do not have this value yet, so store it
        push!(set.vals, val)
        set.slots[ind] = UInt16(length(set.vals))
    end
    set
end

function copyto!(arr::Vector{T}, ts::TupleSet{T}) where {T}
    copyto!(arr, ts.vals)
end

function empty!(ts::TupleSet)
    for i in eachindex(ts.slots)
        ts.slots[i] = zero(UInt16)
    end
    empty!(ts.vals)
    ts
end

function reset!(ts::TupleSet)
    for i in eachindex(ts.slots)
        ts.slots[i] = zero(UInt16)
    end
    ts
end