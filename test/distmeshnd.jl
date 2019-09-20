dim=3
d(p) = sqrt(sum(p.^2,2))-1
p,t = distmeshnd(d,huniform,0.2,[-ones(1,dim);ones(1,dim)],[])