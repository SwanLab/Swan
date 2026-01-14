function [shapef,dershapef,COOR] = QuadrilateralGen_givenPoints(nnodeE,xi) ;
%  Shape functions and derivatives of shape functions for any quadrilateral
%  element 
% Joaquin A. Hernandez,  10th-April-2020.  (28th day of COVID-19 quarantine )
%-----------------------------------------------------------------------------

%[N,B] = Hexahedra27points_symbolic ;

%addpath(genpath('/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/SVDlibrary')) ;


ndim = 2;
% COORDINATES OF THE NODES
% ------------------------
% axy = 1;
% [COOR]= COOR_Quad_n9(axy) ; 
xlim = [-1,1] ; 
ylim = [-1,+1] ; 
nnode_1dir =  sqrt(nnodeE) ; 
if ceil(sqrt(nnode_1dir)) ~= sqrt(nnode_1dir)
    error('Wrong number of nodes') 
end

x = linspace(xlim(1),xlim(2),nnode_1dir) ; 
y = linspace(ylim(1),ylim(2),nnode_1dir) ; 

[x,y] = meshgrid(x,y) ; 

COOR = [x(:),y(:)] ; 


% COORDINATES OF THE  POINT
 
posgp = xi' ; 
% Now we wish to write  the shape function of node I as
% NI = (x^0*y^0*z^0)*a_1 +(x^1*y^0*z^0)*a_2 + .... (x^2*y^2*z^2)*a_27
%
ORDER = nnode_1dir-1; 
n = [ORDER ORDER] ;
[Pnodes ]= CoordinateMatrixPolynomial(COOR,n) ;
COEFFSpol = inv(Pnodes) ;

% Shape functions at the Gauss Points
[Pgauss,Pder ]= CoordinateMatrixPolynomial(posgp,n) ;
shapef = Pgauss*COEFFSpol ;    % ngaus x nnodeE
 ngaus = size(posgp,1) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;

for idim = 1:ndim
   
    dershapef(idim,:,:) =  (Pder{idim}*COEFFSpol)' ;
    
end

 


end

