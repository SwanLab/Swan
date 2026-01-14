function [celasINV3D celas] = Compliance_ElasticityANGLE(E,nu,G,typePROBLEM,ANGLE)

if nargin == 0
    load('tmp1.mat')
end

%%%%
% Compliance matrix for a transverse isotropic material
if length(E)==1
    E = [E E ];
    nu = [nu nu ] ;
    G = [G G ] ;
    
end
%----------------------


%%%%


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
        rowcol = [1 2 6] ;
        celas = celas3D(rowcol,rowcol) ;
    case 'pstress'
        rowcol = [1 2 6] ;
        celasINV = celasINV3D(rowcol,rowcol) ;
        celas = inv(celasINV) ;
    case '3D'
        celas = inv(celasINV3D)  ;
        
        
        if ANGLE ~=0
            T = RotationMatrix(ANGLE) ;
            celas = T*celas*T' ;
            celasINV3D = inv(celas) ;
        end
        
end

