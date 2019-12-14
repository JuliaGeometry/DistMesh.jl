function delaunayn(points)
    # M - no merge of close facets
    # J - no jettison of points
    # B - No boundary info
    # N - No node output
    # F - No face and edge info
    # I - No mesh iteration numbers
    # Q - Quiet
    tetio = tetrahedralize(TetGen.TetgenIO(points), "JBNFIQ")
    tetio
end

function delaunayn_nosort(points)
    tetio = tetrahedralize(TetGen.TetgenIO(points), "Qb/1") # Q- Quiet
    tetio
end

function delaunayn!(fdist, p, t, geps, sorted_pts)
    triangulation = sorted_pts ? delaunayn_nosort(p) : delaunayn(p)
    t_d = triangulation.tetrahedra
    resize!(t, length(t_d))
    copyto!(t, t_d) # we need to copy since we have a shared reference with tetgen

    # average points to get mid point of each tetrahedra
    # if the mid point of the tetrahedra is outside of
    # the boundary we remove it.
    # TODO: this is an inlined filter call. Would be good to revert
    # TODO: can we use the point distance array to pass boundary points to
    #        tetgen so this call is no longer required?
    j = firstindex(t)
    for ai in t
        t[j] = ai
        pm = (p[ai[1]].+p[ai[2]].+p[ai[3]].+p[ai[4]])./4
        j = ifelse(fdist(pm) <= -geps, nextind(t, j), j)
    end
    j <= lastindex(t) && resize!(t, j-1)
    nothing
end

# needs some tweaks, gives garbage results, might need to twek the julia wrapper?
function reconstruct!(points, tets)
    tetio = tetrahedralize(TetGen.TetgenIO(points,tetrahedrons=tets), "Qr") # Q- Quiet, r- retriangulate
    tetio
end
