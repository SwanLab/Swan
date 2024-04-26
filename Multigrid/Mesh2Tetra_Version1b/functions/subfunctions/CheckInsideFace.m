function inter=CheckInsideFace(VN0,VN1,VN2,rx,ry,mi,ma,ignore_corners)
if((rx<mi(1))||(rx>ma(1))||(ry<mi(2))||(ry>ma(2))), inter=false; return, end

Lambda=BarycentricCoordinatesTriangle(VN0,VN1,VN2,rx,ry);

% Check if inside face
if(ignore_corners), mv=1e-8; vv=1-1e-8; else mv=0; vv=1; end
inter=(Lambda(1)>=mv)&&(Lambda(1)<=vv)&&(Lambda(2)>=mv)&&(Lambda(2)<=vv)&& (Lambda(3)>=mv)&&(Lambda(3)<=vv);
