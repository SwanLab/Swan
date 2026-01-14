function [shapef,dershapef,COOR] = Quadrilateral4NInPoints_givenPoints(xi) ;
% This function returns, for each 9-node quadrilateral elements,
% and using a 3x3 Gauss rule (ngaus=0):
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


nnodeE= 4;
ndim = 2;
% COORDINATES OF THE NODES
% ------------------------
COOR = [-1,-1
    +1,-1
    +1,+1
    -1,+1] ;
% COORDINATES OF THE  POINT

posgp = xi' ;
% Now we wish to write  the shape function of node I as
% NI = (x^0*y^0*z^0)*a_1 +(x^1*y^0*z^0)*a_2 + .... (x^2*y^2*z^2)*a_27
%
n = [1 1] ;
[Pnodes ]= CoordinateMatrixPolynomial(COOR,n) ;
COEFFSpol = inv(Pnodes) ;

% Shape functions at the Gauss Points
[Pgauss,Pder ]= CoordinateMatrixPolynomial(posgp,n) ;
shapef = Pgauss*COEFFSpol ;    % ngaus x nnodeE
ndim = 2;  nnodeE = 4 ; ngaus = size(posgp,1) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;

for idim = 1:ndim
    
    dershapef(idim,:,:) =  (Pder{idim}*COEFFSpol)' ;
    
end




end

