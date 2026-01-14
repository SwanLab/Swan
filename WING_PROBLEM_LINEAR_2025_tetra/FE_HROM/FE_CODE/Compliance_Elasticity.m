function [celasINV3D celas] = Compliance_Elasticity(E,nu,G,typePROBLEM,StrainStressWith4Components)

if nargin == 0
    load('tmp1.mat')
elseif nargin == 4 
    StrainStressWith4Components = 0; 
end

%%%%
% Compliance matrix for a transverse isotropic material
if length(E)==1
    E = [E E ];
    nu = [nu nu ] ;
    G = [G G ] ;
    
end
%----------------------
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  G = E/2/(1+nu) ;
celasINV3D = [1/E(1)  -nu(1)/E(1) -nu(1)/E(1)  0 0 0
    -nu(1)/E(1)  1/E(2)  -nu(2)/E(2)  0 0 0
    -nu(1)/E(1)  -nu(2)/E(2)  1/E(2)  0 0 0
    0      0        0  1/G(2)   0 0
    0      0        0  0 1/G(1) 0
    0      0      0  0  0  1/G(1)] ;
switch typePROBLEM
    case 'pstrain'
        celas3D = inv(celasINV3D)  ;
        if StrainStressWith4Components == 1
            rowcol = [1 2  6  3] ;  % Include also the normal stress along the z direction
        else
            rowcol = [1 2 6] ;
        end
        celas = celas3D(rowcol,rowcol) ;
    case 'pstress'
        rowcol = [1 2 6] ;
        celasINV = celasINV3D(rowcol,rowcol) ;
        celas = inv(celasINV) ;
    case '3D'
        celas = inv(celasINV3D)  ;
end