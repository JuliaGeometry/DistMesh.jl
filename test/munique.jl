
a = [1 2; 3 4; 5 6; 1 2; 3 4; 8 8]

@test DistMesh.munique(a) == [1,2,3,1,2,4]