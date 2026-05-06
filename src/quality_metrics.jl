################################################################################
### Element Qualities - Simplex Elements
################################################################################

# Metric 1: Radius Ratio
element_quality_radius_ratio(::Simplex{1}, el) = 1.0

function element_quality_radius_ratio(::Simplex{2}, el)
    p1, p2, p3 = el
    a, b, c = norm(p2 - p1), norm(p3 - p2), norm(p1 - p3)
    s = (a + b + c) / 2
    denom = a * b * c
    return denom ≈ 0 ? 0.0 : 8 * (s - a) * (s - b) * (s - c) / denom
end

function element_quality_radius_ratio(::Simplex{3}, el)
    p1, p2, p3, p4 = el
    d12 = p2 - p1;  d13 = p3 - p1;  d14 = p4 - p1
    d23 = p3 - p2;  d24 = p4 - p2;  d34 = p4 - p3

    v  = element_volume(Simplex{3}(), el)
    s1 = norm(d12 × d13) / 2
    s2 = norm(d12 × d14) / 2
    s3 = norm(d13 × d14) / 2
    s4 = norm(d23 × d24) / 2

    p_1 = norm(d12) * norm(d34)
    p_2 = norm(d23) * norm(d14)
    p_3 = norm(d13) * norm(d24)

    denom = (s1 + s2 + s3 + s4) * sqrt((p_1 + p_2 + p_3) * (p_1 + p_2 - p_3) *
                                        (p_1 + p_3 - p_2) * (p_2 + p_3 - p_1))
    return denom ≈ 0 ? 0.0 : 216 * v^2 / denom
end

# Metric 2: Mean Ratio
element_quality_mean_ratio(::Simplex{1}, el) = 1.0

function element_quality_mean_ratio(::Simplex{2}, el)
    p1, p2, p3 = el
    area = element_volume(Simplex{2}(), el)
    l_sq = sum(abs2, p2 - p1) + sum(abs2, p3 - p2) + sum(abs2, p1 - p3)
    return l_sq ≈ 0 ? 0.0 : (4 * sqrt(3) * area) / l_sq
end

function element_quality_mean_ratio(::Simplex{3}, el)
    p1, p2, p3, p4 = el
    v    = element_volume(Simplex{3}(), el)
    l_sq = sum(abs2, p2 - p1) + sum(abs2, p3 - p1) + sum(abs2, p4 - p1) +
           sum(abs2, p3 - p2) + sum(abs2, p4 - p2) + sum(abs2, p4 - p3)
    return l_sq ≈ 0 ? 0.0 : 216 * v / sqrt(3) / l_sq^(3/2)
end

################################################################################
### Element Qualities - Block Elements
################################################################################

function element_quality_mean_ratio(::Block{2}, el)
    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4)

    l2 = SVector(sum(abs2, e[1]), sum(abs2, e[2]), sum(abs2, e[3]), sum(abs2, e[4]))
    A  = SVector(_cross_mag(e[1], -e[4]), _cross_mag(e[2], -e[1]),
                 _cross_mag(e[3], -e[2]), _cross_mag(e[4], -e[3]))
                 
    Q  = SVector((l2[1]+l2[3]), (l2[2]+l2[4]),
                 (l2[3]+l2[1]), (l2[4]+l2[2])) ./ (2 .* A)

    return 4 / sum(Q)
end

function element_quality_condition_number(::Block{2}, el)
    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4) 

    l2   = SVector(sum(abs2, e[1]), sum(abs2, e[2]), sum(abs2, e[3]), sum(abs2, e[4]))
    sins = SVector(_cross_mag(e[1], -e[4]), _cross_mag(e[2], -e[1]),
                   _cross_mag(e[3], -e[2]), _cross_mag(e[4], -e[3])) ./
           SVector(sqrt(l2[4]*l2[1]), sqrt(l2[1]*l2[2]),
                   sqrt(l2[2]*l2[3]), sqrt(l2[3]*l2[4]))
                   
    any(<=(0), sins) && return 0.0

    k = SVector(l2[4]+l2[1], l2[1]+l2[2], l2[2]+l2[3], l2[3]+l2[4]) ./
        (SVector(sqrt(l2[4]*l2[1]), sqrt(l2[1]*l2[2]),
                 sqrt(l2[2]*l2[3]), sqrt(l2[3]*l2[4])) .* sins)

    return 4 / sqrt(sum(abs2, k))
end

function element_quality_min_scaled_jacobian(::Block{2}, el)
    length(first(el)) == 2 || error("Minimum scaled jacobian for Block{2} is strictly defined for 2D meshes. Surface quad implementation (3D) is not supported.")

    p1, p2, p3, p4 = el
    e = (p2-p1, p3-p2, p4-p3, p1-p4)  
    cross2d(u, v) = u[1]*v[2] - u[2]*v[1]

    J = SVector(cross2d(e[1],-e[4]), cross2d(e[2],-e[1]),
                cross2d(e[3],-e[2]), cross2d(e[4],-e[3]))
    maxJ = maximum(abs.(J))
    return maxJ ≈ 0 ? 0.0 : minimum(J) / maxJ
end

function element_quality_mean_ratio(::Block{3}, el)
    p1,p2,p3,p4,p5,p6,p7,p8 = el
    corners = (
        (p2-p1, p4-p1, p5-p1), (p3-p2, p1-p2, p6-p2),
        (p4-p3, p2-p3, p7-p3), (p1-p4, p3-p4, p8-p4),
        (p8-p5, p6-p5, p1-p5), (p5-p6, p7-p6, p2-p6),
        (p6-p7, p8-p7, p3-p7), (p7-p8, p5-p8, p4-p8),
    )
    function corner_quality(e1, e2, e3)
        detW = dot(e1, cross(e2, e3))
        detW <= 0 && return 0.0
        frob2 = sum(abs2, e1) + sum(abs2, e2) + sum(abs2, e3)
        return 3 * cbrt(detW^2) / frob2
    end
    return minimum(corner_quality(c...) for c in corners)
end

################################################################################
### Default Quality Metrics
################################################################################

default_quality_metric(::Simplex) = element_quality_radius_ratio
default_quality_metric(::Block)   = element_quality_mean_ratio

################################################################################
### User-Facing Shorthands
################################################################################

"""
    element_qualities(m::DMesh; metric=default_quality_metric(E()))

Return a `Vector` of quality metrics for every element in the mesh. 
"""
function element_qualities(m::DMesh{D, T, E}; metric=default_quality_metric(E())) where {D, T, E}
    return collect(element_map(metric, m))
end
