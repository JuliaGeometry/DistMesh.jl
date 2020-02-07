# This file is licensed under the MIT "Expat" License.

# implementing scale-free Hilbert ordering. Real all about it here:
# http://doc.cgal.org/latest/Spatial_sorting/index.html

# original implementation in:
# https://github.com/JuliaGeometry/GeometricalPredicates.jl
# The GeometricalPredicates.jl package is licensed under the MIT "Expat" License:
#    Copyright (c) 2014: Ariel Keselman.

# modifications for StaticArrays: sjkelly (Under terms of MIT License)

const coordinatex = 1
const coordinatey = 2
const coordinatez = 3
next2d(c) = c % 2 + 1
next3d(c) = c % 3 + 1
nextnext3d(c) = (c + 1) % 3 + 1

const forward = true
const backward = false

function select!(direction, coord, v::Array{T,1}, k::Integer, lo::Integer, hi::Integer, carry::CT) where {T<:AbstractVector, CT}
    #lo <= k <= hi || error("select index $k is out of range $lo:$hi")
    if direction == forward
        @inbounds while lo < hi
            if isone(hi-lo)
                if v[hi][coord] < v[lo][coord]
                    v[lo], v[hi] = v[hi], v[lo]
                    if CT !== Nothing; carry[lo], carry[hi] = carry[hi], carry[lo]; end
                end
                return #v[k]
            end
            pivot = v[(lo+hi)>>>1]
            i, j = lo, hi
            while true
                pivot_elt = pivot[coord]
                while v[i][coord] < pivot_elt; i += 1; end
                while pivot_elt < v[j][coord]; j -= 1; end
                i <= j || break
                v[i], v[j] = v[j], v[i]
                if CT !== Nothing; carry[i], carry[j] = carry[j], carry[i]; end
                i += 1; j -= 1
            end
            if k <= j
                hi = j
            elseif i <= k
                lo = i
            else
                return #pivot
            end
        end
    else
        @inbounds while lo < hi
            if isone(hi-lo)
                if v[hi][coord] > v[lo][coord]
                    v[lo], v[hi] = v[hi], v[lo]
                    if CT !== Nothing; carry[lo], carry[hi] = carry[hi], carry[lo]; end
                end
                return #v[k]
            end
            pivot = v[(lo+hi)>>>1]
            i, j = lo, hi
            while true
                pivot_elt = pivot[coord]
                while v[i][coord] > pivot_elt; i += 1; end
                while pivot_elt > v[j][coord]; j -= 1; end
                i <= j || break
                v[i], v[j] = v[j], v[i]
                if CT !== Nothing; carry[i], carry[j] = carry[j], carry[i]; end
                i += 1; j -= 1
            end
            if k <= j
                hi = j
            elseif i <= k
                lo = i
            else
                return #pivot
            end
        end
    end
    #return v[lo]
    nothing
end


# 2D version
# function hilbertsort!(directionx::AbstractDirection, directiony::AbstractDirection, coordinate::AbstractCoordinate, a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64=4) where T<:AbstractPoint2D
#     hi-lo <= lim && return a

#     i2 = (lo+hi)>>>1
#     i1 = (lo+i2)>>>1
#     i3 = (i2+hi)>>>1

#     select!(directionx, coordinate, a, i2, lo, hi)
#     select!(directiony, next2d(coordinate), a, i1, lo, i2)
#     select!(!directiony, next2d(coordinate), a, i3, i2, hi)

#     hilbertsort!(directiony, directionx, next2d(coordinate), a, lo, i1, lim)
#     hilbertsort!(directionx, directiony, coordinate, a, i1, i2, lim)
#     hilbertsort!(directionx, directiony, coordinate, a, i2, i3, lim)
#     hilbertsort!(!directiony, !directionx, next2d(coordinate), a, i3, hi, lim)

#     return a
# end

function hilbertsort!(directionx, directiony, directionz, coordinate, a::Vector, lo::Integer, hi::Integer, lim::Integer, carry)
    hi-lo <= lim && return a

    i4 = (lo+hi)>>>1
    i2 = (lo+i4)>>>1
    i1 = (lo+i2)>>>1
    i3 = (i2+i4)>>>1
    i6 = (i4+hi)>>>1
    i5 = (i4+i6)>>>1
    i7 = (i6+hi)>>>1

    select!(directionx, coordinate, a, i4, lo, hi, carry)
    select!(directiony, next3d(coordinate), a, i2, lo, i4, carry)
    select!(directionz, nextnext3d(coordinate), a, i1, lo, i2, carry)
    select!(!directionz, nextnext3d(coordinate), a, i3, i2, i4, carry)
    select!(!directiony, next3d(coordinate), a, i6, i4, hi, carry)
    select!(directionz, nextnext3d(coordinate), a, i5, i4, i6, carry)
    select!(!directionz, nextnext3d(coordinate), a, i7, i6, hi, carry)

    hilbertsort!( directionz,  directionx,  directiony, nextnext3d(coordinate), a, lo, i1, lim, carry)
    hilbertsort!( directiony,  directionz,  directionx, next3d(coordinate),     a, i1, i2, lim, carry)
    hilbertsort!( directiony,  directionz,  directionx, next3d(coordinate),     a, i2, i3, lim, carry)
    hilbertsort!( directionx, !directiony, !directionz, coordinate,             a, i3, i4, lim, carry)
    hilbertsort!( directionx, !directiony, !directionz, coordinate,             a, i4, i5, lim, carry)
    hilbertsort!(!directiony,  directionz, !directionx, next3d(coordinate),     a, i5, i6, lim, carry)
    hilbertsort!(!directiony,  directionz, !directionx, next3d(coordinate),     a, i6, i7, lim, carry)
    hilbertsort!(!directionz, !directionx,  directiony, nextnext3d(coordinate), a, i7, hi, lim, carry)

    return a
end

#hilbertsort!(a::Array{T,1}) where {T<:AbstractPoint2D} = hilbertsort!(backward, backward, coordinatey, a, 1, length(a))
#hilbertsort!(a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64) where {T<:AbstractPoint2D} = hilbertsort!(backward, backward, coordinatey, a, lo, hi, lim)
"""
Hilbert Sorting. If `carry` is specified, this array will be permuted in line with the
specified array.
"""
hilbertsort!(a::Vector, carry=nothing) = hilbertsort!(backward, backward, backward, coordinatez, a, 1, length(a), 8, carry)
hilbertsort!(a::Vector, lo::Int64, hi::Int64, lim::Int64) = hilbertsort!(backward, backward, backward, coordinatey, a, lo, hi, lim)
