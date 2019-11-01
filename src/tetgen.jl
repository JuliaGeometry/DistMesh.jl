function delaunayn(points)
    tetio = tetrahedralize(TetGen.TetgenIO(points), "Q") # Q- Quiet
    tetio
end

function reconstruct!(points, tets)
    tetio = tetrahedralize(TetGen.TetgenIO(points,tetrahedrons=tets), "Qr") # Q- Quiet, r- retriangulate
    tetio
end
