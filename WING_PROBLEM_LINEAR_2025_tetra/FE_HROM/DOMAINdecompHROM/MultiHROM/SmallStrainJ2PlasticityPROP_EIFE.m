function [MATPRO,DATA] = SmallStrainJ2PlasticityPROP_EIFE(MESH,typePROBLEM,PROPMAT,DATA)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/05_VECTORIZ/02_PLASTICITY.mlx
%
if nargin == 0
    load('tmp.mat')
end

ndim = size(MESH.COOR,2)  ;
if ndim==2
    nstrain = 4;
else
    error('Option not implemented')
end

switch typePROBLEM
    case 'pstrain'
    otherwise
        error('Option not implemented')
end

nelem = size(MESH.MaterialType,1) ;
% MATPRO.celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
MATPRO.dens = zeros(nelem*DATA.MESH.ngaus_RHS,1) ;

%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
    ELEMS = find(MESH.MaterialType == imat) ; % Elements of this type of material
    PROPMATLOC = PROPMAT(imat).PROPMAT ;
    MaterialTypeLocal = PROPMAT(imat).EIFE_prop.BodyForces.MaterialType ;
    dens_elem = AssignOneEIFelementDENSmat_j2(DATA,PROPMATLOC,MaterialTypeLocal)  ;
    ELEMSdofsDENS = small2large(ELEMS,size(dens_elem,1)) ;
    MATPRO.dens(ELEMSdofsDENS,:) = repmat(dens_elem,length(ELEMS),1)  ;
end


DATA = DefaultField(DATA,'StrainStressWith4Components',1) ;
DATA.nstrain = nstrain ;
DATA.ISPSTRAIN = 1;



ngausE = DATA.MESH.ngaus_STRESS ;
DATA.MESH.ngaus = ngausE ; 

ngaus =  nelem*ngausE ;
E_young = zeros(ngaus,1) ;  % Young's modulus
Hmodulus = zeros(ngaus,1) ; % Hardening parameters
PoissonCOEF = zeros(ngaus,1) ;  % Poisson's ratio
sigmay_0 = 1e20*ones(ngaus,1) ;  % Yield stress
%dbstop('30')

for imat = 1:length(PROPMAT)
    %   MATLOC = PROPMAT(imat)  ;
    ELEMS = find(MESH.MaterialType == imat) ;
    PROPMATLOC = PROPMAT(imat).PROPMAT ;
    MaterialTypeLocal = PROPMAT(imat).EIFE_prop.INTforces.MaterialType ;
    
    ELEMSdofs = small2large(ELEMS,ngausE) ;
    
    
     
    
    E_young_loc = AssignOneEIFelement_PROP(ngausE,PROPMATLOC,MaterialTypeLocal,'E') ;
    E_young(ELEMSdofs) = repmat(E_young_loc,length(ELEMS),1)  ;
    
    PoissonCOEF_loc = AssignOneEIFelement_PROP(ngausE,PROPMATLOC,MaterialTypeLocal,'nu') ;
    PoissonCOEF(ELEMSdofs) = repmat(PoissonCOEF_loc,length(ELEMS),1)  ;
    
    sigmay_0_loc = AssignOneEIFelement_PROP(ngausE,PROPMATLOC,MaterialTypeLocal,'sigmay') ;
    sigmay_0(ELEMSdofs) = repmat(sigmay_0_loc,length(ELEMS),1)  ;
    
    Hmodulus_loc = AssignOneEIFelement_PROP(ngausE,PROPMATLOC,MaterialTypeLocal,'Hmodul') ;
    Hmodulus(ELEMSdofs) = repmat(Hmodulus_loc,length(ELEMS),1)  ;
    
    %
    %
    %
    %     switch PROPMAT(imat).TYPECONST
    %         case {'ELASTPLAS','ELAST'}
    %             E_young(pointsMAT) = PROPMAT(imat).E;
    %             PoissonCOEF(pointsMAT) = PROPMAT(imat).nu;
    %             switch  PROPMAT(imat).TYPECONST
    %                 case 'ELASTPLAS'
    %                     switch   PROPMAT(imat).ELASTPLAS.TYPE
    %                         case 'LINEARISOTRO'
    %                             sigmay_0(pointsMAT) =  PROPMAT(imat).sigmay ;
    %                             Hmodulus(pointsMAT) = PROPMAT(imat).Hmodul;
    %                         otherwise
    %                             error('This problem is not amenable to vectorization (stress calculation)')
    %                     end
    %
    %             end
    %     end
end




% Elastic coefficeints (kappa, mu). For all gauss points
ngausT = nstrain*length(E_young) ;
mu=E_young./(2*(1+PoissonCOEF)); % Shear modulus
kappa  = E_young./(3*(1-2*PoissonCOEF)); % Bulk modulus
%muA =reshape([mu'; mu'; mu'; mu'],ngausT,1) ; % For all gauss points and components
%HmodulusA = reshape([Hmodulus'; Hmodulus'; Hmodulus'; Hmodulus'],ngausT,1) ;
%kappaA =reshape([kappa'; kappa'; kappa'; kappa'],ngausT,1) ;
lambda  =  kappa-2*mu/3 ;  % Lame parameter
%lambdaA  =  kappaA-2*muA/3 ;


MATPRO.sigmay_0 = sigmay_0 ;
MATPRO.Hmodulus = Hmodulus ;
MATPRO.mu = mu ;
MATPRO.kappa = kappa ;
%PROPMAT.muA = muA ;
%PROPMAT.kappaA = kappaA ;
%PROPMAT.lambdaA = lambdaA ;
%PROPMAT.HmodulusA = HmodulusA ;



%%%%%%% ELASTIC TANGENT TENSOR (ngaus x 4 matrix)

ident = [1 1 0 1]' ;
volM = ident*ident' ;
identVOL = repmat(volM,ngaus,1) ;
MATPRO.identVOL = identVOL ;

identV = eye(nstrain) ;
identV(3,3) = 0.5 ;
identV = repmat(identV,ngaus,1) ;
MATPRO.identV = identV ;
% ----
% Volumetric term
lambda = kappa-2*mu/3 ;
lambdaA =reshape([lambda'; lambda'; lambda'; lambda'],ngausT,1) ; % For all gauss points and components
muA =reshape([mu'; mu'; mu'; mu'],ngausT,1) ;
C_vol =  bsxfun(@times,identVOL,lambdaA) ;
% Deviatoric term
C_dev =  bsxfun(@times,identV,2*muA) ;


Celas = C_vol+C_dev ;
MATPRO.Celas = Celas ;
