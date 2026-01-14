function [MESH,wSTs,xGAUSS] = GetMeshVariablesGAUSStri(DATA) 

if nargin == 0
    load('tmp.mat')
end

% --------------------------------------
% READING THE MESH 
% -----------------
DATA = DefaultField(DATA,'NameMesh',[]) ; 
if isempty(DATA.NameMesh)
    MESH = DATA.MESH ; 
else
[MESH]= ReadMeshFileStr(DATA.NameMesh,'READ_MATERIAL_COLUMN',1)  ; 
end
% Gauss Points per element 

% if length(DATA.NumberOfGaussPointsPerElement) == 2
%     
%     [~, ~, xGAUSS, weights]  =TensorProd2Ddiscr(DATA.NumberOfGaussPointsPerElement) ;
%     [xx,yy]  = meshgrid(xGAUSS{1},xGAUSS{2}) ;
%     xx = xx(:) ; yy = yy(:);
%     posgp= [xx'; yy'] ;
%     
% elseif length(DATA.NumberOfGaussPointsPerElement) == 3
%     [~, ~, xGAUSS, weights]  =TensorProd3Ddiscr(DATA.NumberOfGaussPointsPerElement) ;
%     [xx,yy,zz]  = meshgrid(xGAUSS{1},xGAUSS{2},xGAUSS{3}) ;
%     xx = xx(:) ; yy = yy(:); zz = zz(:);
%     posgp= [xx'; yy';zz'] ;
% else
%     error('Option not implemented')
%     
% end

% if size(MESH.COOR,1) == 3
%     posgp =  [1/3,1/3] ; 
%     weights = 1/2 ; 
% end
% 
% 
% 
% DATA.posgp_given  = posgp ;
% DATA.weights_given = weights ; 
% % Element shape functions/Derivatives  
% TypeIntegrand = {DATA.posgp_given,DATA.weights_given} ;  
 nnodeE = size(MESH.CN,2) ;

TypeIntegrand = 'K' ; 


[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(MESH.TypeElement,nnodeE,TypeIntegrand) ;
 
[nnode,ndim] = size(MESH.COOR) ; 
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