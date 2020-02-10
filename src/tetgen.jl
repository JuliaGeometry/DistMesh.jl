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

# needs some tweaks, gives garbage results, might need to twek the julia wrapper?
function reconstruct!(points, tets)
    tetio = tetrahedralize(TetGen.TetgenIO(points,tetrahedrons=tets), "Qr") # Q- Quiet, r- retriangulate
    tetio
end
