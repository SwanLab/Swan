function   NormalsPLOT(DATAIN,TypeElementPRINT,NAMEMESH,NameFile_res,NORMALS,TANG2,TANG3)

if nargin ==0
    load('tmp.mat')
end

% Rotation matrix: [NORMALS,TANG2,TANG3]
% -------------------------------------------------------------------------------------------------
ndim = 3; 
NORMALS = rotMATprint(:,1:3:end)' ; 
TANG2 = rotMATprint(:,2:3:end)' ; 
TANG3 = rotMATprint(:,3:3:end)' ; 



  
ElementsPrint =  Indices1DprintRotation' ;
 

%end

nelemE = [] ;
posgp = cell(2,1) ;
NODAL_VECTOR =[] ;
NODAL_SCALAR = [] ;
imesh = 1;
MESH(imesh).NAMEMESH = NAMEMESH{imesh};

MESH(imesh).GAUSS_VECTOR(1).NAME = 'xLOCAL' ;
MESH(imesh).GAUSS_VECTOR(1).COMP = {'x','y','z'}  ;
% All elements but joints
MESH(imesh).GAUSS_VECTOR(1).ELEMENTS =  ElementsPrint';
MESH(imesh).GAUSS_VECTOR(1).VAR = NORMALS' ;

MESH(imesh).GAUSS_VECTOR(2).NAME = 'yLOCAL' ;
MESH(imesh).GAUSS_VECTOR(2).COMP = {'x','y','z'}  ;
% All elements but joints
MESH(imesh).GAUSS_VECTOR(2).ELEMENTS =  ElementsPrint';
MESH(imesh).GAUSS_VECTOR(2).VAR = TANG2' ;

MESH(imesh).GAUSS_VECTOR(3).NAME = 'zLOCAL' ;
MESH(imesh).GAUSS_VECTOR(3).COMP = {'x','y','z'}  ;
% All elements but joints
MESH(imesh).GAUSS_VECTOR(3).ELEMENTS =  ElementsPrint';
MESH(imesh).GAUSS_VECTOR(3).VAR = TANG3' ;

GidResults2DFE_multi(NameFile_res,ndim,NODAL_VECTOR,NODAL_SCALAR,MESH);