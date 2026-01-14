function [N,B] = Hexahedra27points_symbolic 


syms x y z

addpath(genpath('/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/SVDlibrary')) ;


nnodeE= 27; 
ndim = 3;
COOR = zeros(nnodeE,ndim) ; 
% -----------------------------------
% NUMBERING CONVENTION: https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type


% NODES 1-2-3-4: Corners plane z = -1 
z = -1; 
COOR(1,:) = [-1,-1,z] ;
COOR(2,:) = [+1,-1,z] ;
COOR(3,:) = [+1,+1,z] ;
COOR(4,:) = [-1,+1,z] ;
% NODES 5-6-7-8: Corners plane z = +1
z = +1; 
COOR(5,:) = [-1,-1,z] ;
COOR(6,:) = [+1,-1,z] ;
COOR(7,:) = [+1,+1,z] ;
COOR(8,:) = [-1,+1,z] ;
% NODES 9-10-11-12: Midside plane z = -1
z = -1; 
COOR(9,:) = [0,-1,z]  ;
COOR(10,:) = [+1,0,z] ;
COOR(11,:) = [0,+1,z] ;
COOR(12,:) = [-1,0,z] ;
% NODES 13-14-15-16 : Corners plane z = 0 
z = 0; 
COOR(13,:) = [-1,-1,z] ;
COOR(14,:) = [+1,-1,z] ;
COOR(15,:) = [+1,+1,z] ;
COOR(16,:) = [-1,+1,z] ;
% NODES 17-18-19-20: Midside plane z = +1
z =  +1; 
COOR(17,:) = [0,-1,z]  ;
COOR(18,:) = [+1,0,z] ;
COOR(19,:) = [0,+1,z] ;
COOR(20,:) = [-1,0,z] ;
% NODES 21 - Center plane z = -1
z = -1 ; 
COOR(21,:) = [0,0,z] ;
% NODES 22-23-24-25: Midside plane z = +1
z =  0; 
COOR(22,:) = [0,-1,z]  ;
COOR(23,:) = [+1,0,z] ;
COOR(24,:) = [0,+1,z] ;
COOR(25,:) = [-1,0,z] ;
% NODES 26 - Center plane z = +1
z = +1 ; 
COOR(26,:) = [0,0,z] ;
% NODES 27 - Center plane z = 0
z = 0 ; 
COOR(27,:) = [0,0,z] ; 


% Now we wish to write  the shape function of node I as 
% NI = (x^0*y^0*z^0)*a_1 +(x^1*y^0*z^0)*a_2 + .... (x^2*y^2*z^2)*a_27
% 

n = [2 2 2] ; 
[P,Pder ]= CoordinateMatrixPolynomial(COOR,n) ; 

COEFFSpol = inv(P) ; 


% Therefore, each function has 27 coefficients to be fitted. We have 27
% equations. 
 


