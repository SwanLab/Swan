function [MESH,wSTs,xGAUSS] = GetMeshVariables1D(DATA) 

% --------------------------------------
% READING THE MESH 
% -----------------
%[MESH]= ReadMeshFileStr(DATA.NameMesh,'READ_MATERIAL_COLUMN',1)  ; 



MESH.COOR = linspace(DATA.xLIM(1),DATA.xLIM(2),DATA.NumberOfElementsFE+1)' ;
MESH.CN = [(1:DATA.NumberOfElementsFE)',(2:DATA.NumberOfElementsFE+1)'] ; 
MESH.TypeElement = 'Linear'  ; 

 
% Gauss Points per element 
[x1d, weights] = GaussQuad(DATA.NumberOfGaussPointsPerElement, -1, +1) ;


 
posgp= x1d(:)';
DATA.posgp_given  = posgp ;
DATA.weights_given = weights ; 
% Element shape functions/Derivatives  
TypeIntegrand = {DATA.posgp_given,DATA.weights_given} ;  
nnodeE = size(MESH.CN,2) ;  
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(MESH.TypeElement,nnodeE,TypeIntegrand) ;
 
[nnode,ndim] = size(MESH.COOR) ; 
% Next we compute the set of weights and the position of the GAuss points
% for all the elements of the mesh 
[   wSTs,~, posgp_RHS] = ComputeW_RHS(MESH.COOR,MESH.CN,MESH.TypeElement,ndim,TypeIntegrand)  ;




[nelem,nnodeE] = size(MESH.CN) ; 
[ Nshapeelem,~  ] = ComputeNelemALL(MESH.TypeElement,nnodeE,ndim,nelem,TypeIntegrand) ;
ngaus_RHS = size(posgp_RHS,2) ;
MESH.ngausE = ngaus_RHS ; 
MESH.Nst = AssemblyNGlobal(Nshapeelem,nelem,nnodeE,ndim,ngaus_RHS,MESH.CN,nnode) ;
% Integretion points positions 
xNODES= MESH.COOR' ; 
xNODES = xNODES(:) ; 
xGAUSS = MESH.Nst*xNODES; 
 xGAUSS =reshape(xGAUSS,ndim,[])' ; 

% CARTESIAN DOMAIN IN WHICH THE DOMAIN IS EMBEDDED 
MESH.xLIM = zeros(ndim,2) ; 
for idim = 1:ndim 
    MESH.xLIM(idim,1) = min(MESH.COOR(:,idim)) ; 
    MESH.xLIM(idim,2) = max(MESH.COOR(:,idim)) ; 
end