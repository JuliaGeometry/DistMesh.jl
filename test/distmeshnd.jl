dim=3
d(p) = sqrt(sum(p.^2))-1
p,t = distmeshnd(d,huniform,0.2)


#fi = vtk_grid("tet",[x[1] for x in p],[x[2] for x in p],[x[3] for x in p],map(x->MeshCell(10, [x...]),t))
#vtk_cell_data(fi, )
#vtk_save(fi)