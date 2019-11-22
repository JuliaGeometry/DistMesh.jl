function centraldiff(f::Function,p::VT) where VT
    deps = sqrt(eps(eltype(VT)))
    dx = (f(p.+VT(deps,0,0)) - f(p.-VT(deps,0,0)))
    dy = (f(p.+VT(0,deps,0)) - f(p.-VT(0,deps,0)))
    dz = (f(p.+VT(0,0,deps)) - f(p.-VT(0,0,deps)))
    grad = VT(dx,dy,dz)./(2deps) #normalize?
end
