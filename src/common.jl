
function compute_displacements!(fh, dp, pair, L, L0, bars, p, setup,
                                ::Type{VertType}) where {VertType}

    non_uniform = isa(typeof(L0), AbstractVector)

    # compute edge lengths (L) and adaptive edge lengths (L0)
    # Lp norm (p=3) is partially computed here
    Lsum = zero(eltype(L))
    L0sum = non_uniform ? zero(eltype(L0)) : length(pair)
    for i in eachindex(pair)
        pb = pair[i]
        b1 = p[pb[1]]
        b2 = p[pb[2]]
        barvec = b1 - b2 # bar vector
        bars[i] = barvec
        L[i] = sqrt(sum(barvec.^2)) # length
        non_uniform && (L0[i] = fh((b1+b2)./2))
        Lsum = Lsum + L[i].^3
        non_uniform && (L0sum = L0sum + L0[i].^3)
    end

    # zero out force at each node
    for i in eachindex(dp)
        dp[i] = zero(VertType)
    end

    # this is not hoisted correctly in the loop so we initialize here
    # finish computing the Lp norm (p=3)
    lscbrt = (1+(0.4/2^2))*cbrt(Lsum/L0sum)

    # Move mesh points based on edge lengths L and forces F
    for i in eachindex(pair)
        if non_uniform && L[i] < L0[i]*lscbrt || L[i] < lscbrt
            L0_f = non_uniform ? L0[i].*lscbrt : lscbrt
            # compute force vectors
            F = setup.nonlinear ? (L[i]+L0_f)*(L0_f-L[i])/(2*L0_f) : L0_f-L[i]
            # edges are not allowed to pull, only repel
            FBar = bars[i].*F./L[i]
            # add the force vector to the node
            b1 = pair[i][1]
            b2 = pair[i][2]
            dp[b1] = dp[b1] .+ FBar
            dp[b2] = dp[b2] .- FBar
        end
    end
end
