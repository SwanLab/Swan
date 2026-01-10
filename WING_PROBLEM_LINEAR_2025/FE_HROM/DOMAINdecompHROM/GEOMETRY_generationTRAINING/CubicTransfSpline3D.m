function [COORtransf,CentrF1,CentrF2] = ...
    CubicTransfSpline3D(P1_G,P2_G,R_G1,R_12,COORref,COORrefDUMMY,ax,DATAcubic)
% Given the coordinates of the reference domain COORref,
% CubicTransfSpline3D maps this coordinates into a deformed domain
% fulfilling the following conditions 
% 1) Centroid of face 1 coincides with P1_G
% 2) Centrod of face 2 coincides with P2_G
% 3) Midline is normal to the plane at Pini whose normal is R_G1(:,1)
% 4) Midline is normal to the plane at Pfin whose normal is R_12(:,2) =
% (ROTini*dROT)(:,2)
% JAHO, 22-Oct-2020
if nargin == 0
    load('tmp2.mat')
end
dP_G = P2_G-P1_G ;  % Vector joining Pini and Pfin (global coordinates)
% Rotation
dP_1 = R_G1'*dP_G ;   % Vector joining the centroids of the slice, expressed in the coordinates intrinsic to 
                          % centroid 1 
% --------------------------
Dx =  dP_1(1) ; 
Dy = dP_1(2) ; 
Dz = dP_1(3) ; 
t2_1  = R_12(:,1) ;   % Vector tangent to the midline at point 2, expressed in the RS1
% --------------------------------------------------------------------------------------
% System of equations for determining the coefficients of the
% polynomial
% ---------------------------------------------------------------
T = @(lambda) [1,lambda,lambda.^2,lambda.^3]; 
dT =  @(lambda)([0,1,2*lambda,3*lambda.^2]) ; 
COEFF = zeros(8,8) ; 
COEFF(1,1:4) = T(0) ; 
COEFF(2,5:end) = T(0) ; 
COEFF(3,1:4) = dT(0) ; 
COEFF(4,5:end) = dT(0) ; 
COEFF(5,1:4) = T(Dx) ; 
COEFF(6,5:end) = T(Dx) ; 
COEFF(7,1:4) = dT(Dx) ; 
COEFF(8,5:end) = dT(Dx) ; 
b  = zeros(8,1) ; 
b(5) = Dy ; 
b(6) = Dz ; 
b(7) = t2_1(2)/t2_1(1) ; 
b(8) = t2_1(3)/t2_1(1) ; 

a = COEFF\b;
ay = a(1:4) ;
az  = a(5:end) ;
ay = ay([4,3,2,1]) ;  % So that the polynomial  can be evaluated by polyval
az = az([4,3,2,1]) ;

%yini = f(xini) ; 

COORtransfdummy = TransformToCubicRefDomain3D(COORrefDUMMY,dP_1,R_G1,ay,az,P1_G,ax,DATAcubic)   ; 


 

[COORtransf,CentrF1,CentrF2 ]= TransformToCubicRefDomain3D(COORref,dP_1,R_G1,ay,az,P1_G,ax,DATAcubic)   ; 



plot3(COORtransfdummy(:,1),COORtransfdummy(:,2),COORtransfdummy(:,3),'go')