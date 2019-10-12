function delaunayn(points::Vector{Point{3, T}}) where T <: AbstractFloat
    tetio = tetrahedralize(TetGen.TetgenIO(points), "Qe") # Q- Quiet, e - edges
    tetio
end