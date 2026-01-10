function [MATPRO] = NeoHookPROP(MESH,typePROBLEM,PROPMAT,DATA)

if nargin == 0
    load('tmp2.mat')
end

ndim = size(MESH.COOR,2)  ;
if ndim==2
    nstrain = 3;
else
    nstrain = 6 ;
    typePROBLEM ='3D' ;
end

nstrain = DATA.MESH.nstrain  ; 
ngaus = DATA.MESH.ngaus  ; 

nelem = size(MESH.MaterialType,1) ;
%MATPRO.celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
MATPRO.dens = zeros(nelem,1) ;
MATPRO.mu_0 = zeros(nelem*ngaus,1) ;
MATPRO.lambda_0 = zeros(nelem*ngaus,1)  ;



%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
%    celas3D =PROPMAT(imat).ElasticityMatrix ; %
    E = PROPMAT(imat).YoungModulus ; %
    nu = PROPMAT(imat).PoissonRatio ; % 
    
    lambda_0 =  E*nu/(1+nu)/(1-2*nu) ; 
    mu_0 = E/2/(1+nu) ; 
     
 %   INVcelas3D = inv(celas3D) ;
    ELEMS = find(MESH.MaterialType == imat) ;
    
    GAUSS = small2large(ELEMS,ngaus) ; 
    
    MATPRO.mu_0(GAUSS) = mu_0 ; 
    MATPRO.lambda_0(GAUSS) = lambda_0 ; 
    
%     switch typePROBLEM
%         case 'pstrain'
%             rowcol = [1 2 6] ;
%             celas = celas3D(rowcol,rowcol) ;
%         case 'pstress'
%             error('Option not implented yet')
% %             rowcol = [1 2 6] ;
% %             celasINV3D = inv(celas3D) ;
% %             celasINV = celasINV3D(rowcol,rowcol) ;
% %             celas = inv(celasINV) ;
%         case '3D'
%             celas = celas3D ;
%     end
  %  for eLOC=1:length(ELEMS)
       % e = ELEMS(eLOC) ;
     %   MATPRO.celasglo(:,:,e) = celas ;
        MATPRO.dens(ELEMS) = PROPMAT(imat).Density ;
        
   % end
end
