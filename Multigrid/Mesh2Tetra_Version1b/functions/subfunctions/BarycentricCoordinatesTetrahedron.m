function Lambda=BarycentricCoordinatesTetrahedron(p1,p2,p3,p4,rx,ry,rz)
% Edge Vectors
r1 = (p2-p1); r2 = (p3-p1);
r3 = (p4-p1); r4 = (p3-p2);
r5 = (p2-p4); r6 = (p4-p3);

% Jacobian determinant
J = dot(cross(r1,r2),r3); 

% Lambda Gradient 
g1 = cross(r4,r5)/J; 
g2 = cross(r2,r3)/J;
g3 = cross(r3,r1)/J; 
g4 = cross(r1,r2)/J;

% Center compensation
r_b = (p1 + p2 + p3 + p4)/4;
c1=dot(-r_b, g1)+1/4;
c2=dot(-r_b, g2)+1/4;
c3=dot(-r_b, g3)+1/4;
c4=dot(-r_b, g4)+1/4;

% Current location
r=[rx ry rz];
            
% Interpolation values
Lambda = [dot(r,g1)+c1 dot(r,g2)+c2 dot(r,g3)+c3 dot(r,g4)+c4]; 
            