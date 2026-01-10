function [DATA,PROPMAT,densGLO,CelasGLO] = J2initialVARIABLES(nstrain,DATA,nelem,MaterialType,MATERIAL,ngausE)


DATA = DefaultField(DATA,'StrainStressWith4Components',1) ;  

DATA.nstrain = nstrain ; 
EXIST_DENS = 0 ;

DATA.ISPSTRAIN = 0; ;

%celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
%celasgloINV = zeros(6,6,nelem) ;  % Global array of compliance matrices (3D)
densGLO = ones(nelem,1) ;  % Global array of compliance matrices (3D)


NMAT_fe = length(unique(MaterialType));
nmat = length(MATERIAL.PLY) ;

if NMAT_fe > nmat
    for imat = nmat:NMAT_fe
        MATERIAL.PLY(imat) =  MATERIAL.PLY(nmat)   ;
    end
    nmat = NMAT_fe ;
end


if isempty(MaterialType)
    MaterialType = ones(nelem,1) ;
end


ISPSTRAIN = 0 ;

ngaus =  nelem*ngausE ;
E_young = zeros(ngaus,1) ;
Hmodulus = zeros(ngaus,1) ;
PoissonCOEF = zeros(ngaus,1) ;
sigmay_0 = 1e20*ones(ngaus,1) ;
%dbstop('30')
for imat = 1:nmat
    % Load elasticity tensor
    MATLOC = MATERIAL.PLY(imat)  ;
    E = MATLOC.E ;
    nu = MATLOC.nu ;
    G = MATLOC.G ;
    typePROBLEM = MATLOC.typePROBLEM  ;
    StrainStressWith4Components = 0 ;
    switch typePROBLEM
        case '3D'
            error('Option not implemented')
            %                 if ndim ==2
            %                     error('Set typePROBLEM to either pstress or pstrain')
            %                 end
            
        case 'pstrain'
            StrainStressWith4Components = DATA.StrainStressWith4Components ;
            DATA.ISPSTRAIN = 1;
        case 'pstress'
            error('Option not implemented')
            if DATA.StrainStressWith4Components == 1
                error('Incompatible option. ')
            end
    end
    
    DATA.StrainStressWith4Components = StrainStressWith4Components;
    if isfield(MATLOC,'densCOMP')
        dens = MATLOC.densCOMP ; % Density
        EXIST_DENS  =1;
    end
    
    
    elemMAT = find(MaterialType == imat) ;
    pointsMAT = small2large(elemMAT,ngausE) ; 
    
    
    switch MATERIAL.PLY(imat).TYPECONST
        case {'ELASTPLAS','ELAST'}
            E_young(pointsMAT) = MATERIAL.PLY(imat).E;
            PoissonCOEF(pointsMAT) = MATERIAL.PLY(imat).nu;
            switch  MATERIAL.PLY(imat).TYPECONST
                case 'ELASTPLAS'
                    switch   MATERIAL.PLY(imat).ELASTPLAS.TYPE
                        case 'LINEARISOTRO'
                            sigmay_0(pointsMAT) =  MATERIAL.PLY(imat).sigmay ;
                            Hmodulus(pointsMAT) = MATERIAL.PLY(imat).Hmodul;
                            
                        otherwise
                            error('This problem is not amenable to vectorization (stress calculation)')
                    end
                    %   otherwise
                    %      error('This problem is not amenable to vectorization (stress calculation)')
            end
    end
    
    
    
    
    if EXIST_DENS==1
        densGLO(elemMAT) = dens ;
    end
    
    
    
    
    
end




% Elastic coefficeints (kappa, mu). For all gauss points
ngausT = nstrain*length(E_young) ;
mu=E_young./(2*(1+PoissonCOEF));
kappa  = E_young./(3*(1-2*PoissonCOEF));
%muA =reshape([mu'; mu'; mu'; mu'],ngausT,1) ; % For all gauss points and components
%HmodulusA = reshape([Hmodulus'; Hmodulus'; Hmodulus'; Hmodulus'],ngausT,1) ;
%kappaA =reshape([kappa'; kappa'; kappa'; kappa'],ngausT,1) ;
lambda  =  kappa-2*mu/3 ;
%lambdaA  =  kappaA-2*muA/3 ;


PROPMAT.sigmay_0 = sigmay_0 ;
PROPMAT.Hmodulus = Hmodulus ;
PROPMAT.mu = mu ;
PROPMAT.kappa = kappa ;
%PROPMAT.muA = muA ;
%PROPMAT.kappaA = kappaA ;
%PROPMAT.lambdaA = lambdaA ;
%PROPMAT.HmodulusA = HmodulusA ;



%%%%%%% ELASTIC TANGENT TENSOR (ngaus x 4 matrix)
 
ident = [1 1 0 1]' ;
volM = ident*ident' ;
identVOL = repmat(volM,ngaus,1) ;
PROPMAT.identVOL = identVOL ;

identV = eye(nstrain) ;
identV(3,3) = 0.5 ;
identV = repmat(identV,ngaus,1) ;
PROPMAT.identV = identV ;
% ----
% Volumetric term
lambda = kappa-2*mu/3 ;
lambdaA =reshape([lambda'; lambda'; lambda'; lambda'],ngausT,1) ; % For all gauss points and components
muA =reshape([mu'; mu'; mu'; mu'],ngausT,1) ;
C_vol =  bsxfun(@times,identVOL,lambdaA) ;
% Deviatoric term
C_dev =  bsxfun(@times,identV,2*muA) ;


Celas = C_vol+C_dev ;
PROPMAT.Celas = Celas ;

%%% Indices for obtaining Ctang
[indI,indJ] = IndicesCtang(size(Celas,1),size(Celas,2)) ;
PROPMAT.indI = indI ;
PROPMAT.indJ = indJ ;

% Diagonal matrix (Elastic range)
  CelasGLO =     ConvertCmatSparseMatrix(Celas,indI,indJ) ;