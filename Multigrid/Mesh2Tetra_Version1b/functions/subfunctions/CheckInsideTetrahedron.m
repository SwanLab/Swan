function inter=CheckInsideTetrahedron(p1,p2,p3,p4,rx,ry,rz,mi,ma)
if((rx<mi(1))||(rx>ma(1))||(ry<mi(2))||(ry>ma(2))||(rz<mi(3))||(rz>ma(3))), inter=false; return, end

Lambda=BarycentricCoordinatesTetrahedron(p1,p2,p3,p4,rx,ry,rz);

% Check if inside tetrahedron
inter=(Lambda(1)>=0)&&(Lambda(1)<=1)&&(Lambda(2)>=0)&&(Lambda(2)<=1)&&(Lambda(3)>=0)&&(Lambda(3)<=1)&&(Lambda(4)>=0)&&(Lambda(4)<=1);
