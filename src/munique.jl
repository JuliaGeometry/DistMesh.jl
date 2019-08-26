
# function to generate value to index mapping
# not sure if tuple is optimal here, but want to ensure
# immutability of keys
rowpair(a, i) = (a[i,1], a[i,2]) => i 

"""
Compatibility with Matlab-style `[C, ia, ic] = unique(A, 'rows')`.

Defintely not a generic implementation, or optimal
"""
function munique(a)
    ua = unique(a, dims=1)
    @assert size(a)[2] == 2

    d = Dict(rowpair(ua, i) for i=1:size(ua, 1))
    # prealloc
    revind = Vector{Int}(undef, size(a, 1))
    # generate reverse mapping by lookup
    for i = 1:size(a, 1)
        revind[i] = d[(a[i,1], a[i,2])]
    end

    revind
end