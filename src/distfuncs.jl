################################################################################
### 1. Basic Shapes (Circles, Rectangles)
################################################################################

"""
    dhypersphere(p, c, r)

Signed distance function for a hypersphere of radius `r` centered at `c`.
Works for any dimension, provided `p` and `c` have matching lengths.
"""
dhypersphere(p, c, r) = norm(p - c) - r

"""
    dcircle(p; c=(0, 0), r=1.0)

Signed distance function for a 2D circle. 
Wrapper around `dhypersphere` that enforces 2D inputs.
"""
dcircle(p; c=(0, 0), r=1) = dhypersphere(SVector{2}(p), SVector{2}(c), r)

"""
    dsphere(p; c=(0, 0, 0), r=1.0)

Signed distance function for a 3D sphere.
Wrapper around `dhypersphere` that enforces 3D inputs.
"""
dsphere(p; c=(0, 0, 0), r=1) = dhypersphere(SVector{3}(p), SVector{3}(c), r)

"""
    drectangle(p, x1, x2, y1, y2)

Signed distance function for a 2D axis-aligned rectangle.
Returns negative values inside the region `[x1, x2] × [y1, y2]`.
"""
drectangle(p, x1, x2, y1, y2) = -minimum((-y1 + p[2], y2 - p[2], -x1 + p[1], x2 - p[1]))

"""
    dblock(p, x1, x2, y1, y2, z1, z2)

Signed distance function for a 3D axis-aligned block (cuboid).
Returns negative values inside the region `[x1, x2] × [y1, y2] × [z1, z2]`.
"""
dblock(p, x1, x2, y1, y2, z1, z2) = -minimum((-z1 + p[3], z2 - p[3],
                                              -y1 + p[2], y2 - p[2],
                                              -x1 + p[1], x2 - p[1]))

################################################################################
### 2. Polygons & Lines
################################################################################

"""
    dline(p, p1, p2)

Distance from point `p` to the finite line segment defined by endpoints `p1` and `p2`.
Returns the Euclidean distance (always non-negative).
"""
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

"""
    inpolygon(p, pv)

Point-in-polygon test using the Ray Casting algorithm.
Returns `true` if `p` is inside the polygon defined by vertices `pv`.
"""
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

# Internal helper (not typically exported, so no docstring needed unless you want one)
dsegment(p, pv) = (dline(p, pv[i+1], pv[i]) for i in 1:length(pv)-1)

"""
    dpoly(p, pv)

Signed distance function for a polygon defined by vertices `pv`.
**Note:** `pv` must be a closed loop (i.e., `pv[end] == pv[1]`).
Returns negative values inside, positive outside.
"""
function dpoly(p, pv)
    d = minimum(dsegment(p, pv))
    inpolygon(p, pv) ? -d : d
end

################################################################################
### 3. Boolean Operations (CSG)
################################################################################

"""
    ddiff(d1, d2)

Difference of two regions (Region 1 minus Region 2).
`d1` and `d2` are signed distances.
"""
ddiff(d1, d2) = max(d1, -d2)

"""
    dunion(d1, d2)

Union of two regions.
`d1` and `d2` are signed distances.
"""
dunion(d1, d2) = min(d1, d2)

"""
    dintersect(d1, d2)

Intersection of two regions.
`d1` and `d2` are signed distances.
"""
dintersect(d1, d2) = max(d1, d2)

################################################################################
### 4. Special Functions
################################################################################

const naca_coeffs = 0.12 / 0.2 * SVector(0.2969, -0.1260, -0.3516, 0.2843, -0.1036)

"""
    dnaca(p)

Implicit level-set function for a NACA 0012 airfoil.
Zero contour defines the airfoil boundary.
"""
function dnaca(p)
    a0, a14... = naca_coeffs
    (abs(p[2]) - sum(a14 .* p[1] .^ SVector(1, 2, 3, 4))) .^ 2 - a0^2 * p[1]
end
