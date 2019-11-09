function huniform(::T) where {T<:AbstractVector}
    one(eltype(T))
end

"""
    HUniform

Empty struct to allow eliding of adaptive meshing function calls.
Use this in place of an adaptive edge size function for better performance.
"""
struct HUniform end