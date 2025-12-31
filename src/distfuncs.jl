#using StaticArrays
#using LinearAlgebra

dhypersphere(p, c, r) = norm(p - c) - r
dcircle(p; c=(0, 0), r=1) = dhypersphere(SVector{2}(p), SVector{2}(c), r)
dsphere(p; c=(0, 0, 0), r=1) = dhypersphere(SVector{3}(p), SVector{3}(c), r)

drectangle(p, x1, x2, y1, y2) = -minimum((-y1 + p[2], y2 - p[2], -x1 + p[1], x2 - p[1]))
dblock(p, x1, x2, y1, y2, z1, z2) = -minimum((-z1 + p[3], z2 - p[3],
                                              -y1 + p[2], y2 - p[2],
                                              -x1 + p[1], x2 - p[1]))

function dline(p, p1, p2)
    p,p1,p2 = SVector{2}.((p,p1,p2))

    v = p2 - p1
    w = p - p1
    c1 = dot(v, w)
    c2 = dot(v, v)
    if c1 <= 0
        return norm(p - p1)
    elseif c1 >= c2
        return norm(p - p2)
    else
        return norm(p - (p1 + c1 / c2 * v))
    end
end

function inpolygon(p, pv)
    cn = 0 # Crossing number counter
    for i = 1:length(pv)-1 # Loop over all edges of the polygon
        if pv[i][2] <= p[2] && pv[i+1][2] > p[2] ||  # Upward crossing
           pv[i][2] > p[2] && pv[i+1][2] <= p[2]     # Downward crossing
            vt = (p[2] - pv[i][2]) / (pv[i+1][2] - pv[i][2]) # Intersection
            if p[1] < pv[i][1] + vt * (pv[i+1][1] - pv[i][1])
                cn += 1 # A valid cross right of p[1]
            end
        end
    end
    return cn % 2 == 1
end

dsegment(p, pv) = (dline(p, pv[i+1], pv[i]) for i in 1:length(pv)-1)
function dpoly(p, pv)
    d = minimum(dsegment(p, pv))
    inpolygon(p, pv) ? -d : d
end

ddiff(d1, d2) = max(d1, -d2)
dunion(d1, d2) = min(d1, d2)
dintersect(d1, d2) = max(d1, d2)

huniform(p) = 1

const naca_coeffs = 0.12 / 0.2 * SVector(0.2969, -0.1260, -0.3516, 0.2843, -0.1036)

function dnaca(p)
    a0, a14... = naca_coeffs
    (abs(p[2]) - sum(a14 .* p[1] .^ SVector(1, 2, 3, 4))) .^ 2 - a0^2 * p[1]
end
