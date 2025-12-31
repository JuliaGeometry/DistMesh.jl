#import Plots

#=
function simpplot(p, t, tdata=nothing; args...)
    pxy(d) = [ tix==0 ? NaN : p[tt[tix]][d] for tix in (1,2,3,1,0), tt in t ]
    args = (args..., seriestype=:shape, aspect_ratio=:equal, leg=false)
    if isnothing(tdata)
        args = (args..., fillcolor=Plots.RGB(0.8, 0.9, 1.0))
    else
        # todo
        #data = vcat(tdata[t], fill(NaN, 1, size(t, 2)))[:]
        #args = (args..., fillcolor=:viridis, fill_z=data, colorbar=:right)
    end
    h = Plots.plot(vec(pxy(1)), vec(pxy(2)); args...)
end
=#

function simplex_area(el)
    if length(el) == 3 # Triangle
        p12 = el[2] - el[1]
        p13 = el[3] - el[1]
        return (p12[1] * p13[2] - p12[2] * p13[1]) / 2
    else
        @assert "Dimension not implemented"
    end
end


function simplex_qual(el)
    if length(el) == 3 # Triangle
        norm(vec) = sqrt(sum(vec.^2))
        a,b,c = ( norm(el[ix[2]] - el[ix[1]]) for ix in ((1,2),(2,3),(3,1)) )
        r = 0.5*sqrt((b+c-a) * (c+a-b) * (a+b-c) / (a+b+c))
        R = a*b*c / sqrt((a+b+c) * (b+c-a) * (c+a-b) * (a+b-c))
        return 2*r/R
    else
        @assert "Dimension not implemented"
    end
end

elemwise_feval(p,t,f) = [ f(p[tt]) for tt in t ]
simpqual(p, t, fqual=simplex_qual) = elemwise_feval(p, t, fqual)
simpvol(p,t) = elemwise_feval(p, t, simplex_area)

snap(x::T, scaling=1) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling*tol*round(x/scaling/tol) + zero(T)  # Adding zero to uniquify -0.0 and 0.0

function fixmesh(p, t; output_ix=false)
    scaling = maximum(norm.(p))
    pp = [ snap.(p1,scaling) for p1 in p ]
    ppp = unique(pp)
    ix = Int.(indexin(ppp, pp))
    jx = Int.(indexin(pp, ppp))
    p = p[ix]
    t = [jx[tt] for tt in t]
    return output_ix ? (p,t,ix) : (p,t)
end 

