function [shapef,dershapef,COOR] = Hexahedra27NInPoints_givenpoint(xi) ;
% This function returns, for each 27-node quadrilateral elements,
% and using a 3x3x3 Gauss rule (ngaus=0):
% NUMBERING CONVENTION: https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type
% weig = Vector of Gauss weights (1xngaus)
% posgp: Position of Gauss points  (ndim x ngaus)
% shapef: Array of shape functions (ngaus x nnodeE)
% dershape: Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)
% Joaquin A. Hernandez,  10th-April-2020.  (28th day of COVID-19 quarantine )
%-----------------------------------------------------------------------------

%[N,B] = Hexahedra27points_symbolic ;

%addpath(genpath('/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/SVDlibrary')) ;

nnodeE= 27;
ndim = 3;
% COORDINATES OF THE NODES
% ------------------------
axy = 1;
bz = 1;
COOR = COOR_Hexag27(axy,bz) ;
% COORDINATES OF THE  POINT
axy = sqrt(3/5) ;
bz = axy ;
posgp = xi'  ; 

% Now we wish to write  the shape function of node I as
% NI = (x^0*y^0*z^0)*a_1 +(x^1*y^0*z^0)*a_2 + .... (x^2*y^2*z^2)*a_27
%
n = [2 2 2] ;
[Pnodes ]= CoordinateMatrixPolynomial(COOR,n) ;
COEFFSpol = inv(Pnodes) ;

% Shape functions at the Gauss Points
[Pgauss,Pder ]= CoordinateMatrixPolynomial(posgp,n) ;
shapef = Pgauss*COEFFSpol ;    % ngaus x nnodeE
ndim = 3;  nnodeE = 27 ; ngaus = size(posgp,1) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;

for idim = 1:ndim
   
    dershapef(idim,:,:) =  (Pder{idim}*COEFFSpol)' ;
    
end

 


end

