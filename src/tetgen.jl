function delaunayn(points)
    tetio = tetrahedralize(TetGen.TetgenIO(points), "Q") # Q- Quiet
    tetio
end

function delaunaynsort(points)
    tetio = tetrahedralize(TetGen.TetgenIO(points), "Qb") # Q- Quiet, b- BRIO-hilbert sort
    tetio
end

# needs some tweaks, gives garbage results, might need to twek the julia wrapper?
function reconstruct!(points, tets)
    tetio = tetrahedralize(TetGen.TetgenIO(points,tetrahedrons=tets), "Qr") # Q- Quiet, r- retriangulate
    tetio
end
