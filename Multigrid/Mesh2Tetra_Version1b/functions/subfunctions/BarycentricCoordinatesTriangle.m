function  Lambda=BarycentricCoordinatesTriangle(VN0,VN1,VN2,rx,ry)

% Normalization factors
f12 = ( VN1(2) - VN2(2) ) * VN0(1)  + (VN2(1) - VN1(1) ) * VN0(2) + VN1(1) * VN2(2) - VN2(1) *VN1(2);
f20 = ( VN2(2) - VN0(2) ) * VN1(1)  + (VN0(1) - VN2(1) ) * VN1(2) + VN2(1) * VN0(2) - VN0(1) *VN2(2);
f01 = ( VN0(2) - VN1(2) ) * VN2(1)  + (VN1(1) - VN0(1) ) * VN2(2) + VN0(1) * VN1(2) - VN1(1) *VN0(2);

% Lambda Gradient
g12x = ( VN1(2) - VN2(2) )/f12; g12y = ( VN2(1) - VN1(1) )/f12;
g20x = ( VN2(2) - VN0(2) )/f20; g20y = ( VN0(1) - VN2(1) )/f20;
g01x = ( VN0(2) - VN1(2) )/f01; g01y = ( VN1(1) - VN0(1) )/f01;

% Center compensation
c12 = (VN1(1) * VN2(2) - VN2(1) * VN1(2))/f12;
c20 = (VN2(1) * VN0(2) - VN0(1) * VN2(2))/f20;
c01 = (VN0(1) * VN1(2) - VN1(1) * VN0(2))/f01;

% Interpolation values
Lambda=[g12x*rx+g12y*ry+c12 g20x*rx+g20y*ry+c20 g01x*rx+g01y*ry+c01];