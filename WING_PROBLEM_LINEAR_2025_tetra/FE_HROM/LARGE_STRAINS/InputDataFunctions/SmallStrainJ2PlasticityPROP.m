function [MATPRO,DATA] = SmallStrainJ2PlasticityPROP(MESH,typePROBLEM,PROPMAT,DATA)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
 % % =========================================================================
% SmallStrainJ2PlasticityPROP — Build per-Gauss-point elastic props (2D p-strain)
% =========================================================================
% PURPOSE
%   Prepare material property arrays and constant elastic tensors at Gauss
%   points for small-strain 2D problems, currently **plane strain** only.
%   This routine gathers FE-mesh material assignments, expands per-material
%   properties to element/Gauss-point fields, and assembles:
%     - Elastic moduli (E, ν) → (μ, κ, λ)
%     - Initial yield stress and isotropic hardening modulus (for J2 plasticity)
%     - The 4×4 elastic stiffness matrix Celas in Voigt form (per Gauss point)
%
% SCOPE / STATUS
%   - Supported kinematics: 2D **plane strain** with Voigt order [xx, yy, xy, zz]
%     (nstrain = 4, where γ_xy is engineering shear; C(3,3) uses 2μ with 0.5
%     scaling embedded in identV).
%   - Supported constitutive types (per PROPMAT(imat).TYPECONST):
%       'ELAST'      → linear elastic
%       'ELASTPLAS'  → J2 elastoplastic (this routine only sets ELASTIC data
%                      + (σ_y0, H) fields; plastic update is done elsewhere).
%   - Multiple materials handled via MESH.MaterialType (int labels per element).
%
% INPUTS
%   MESH : struct with
%            .COOR          → (nnode × ndim) nodal coordinates
%            .MaterialType  → (nelem × 1) integer material label per element
%            .ngaus         → Gauss points per element (DATA.MESH.ngaus is used)
%   typePROBLEM : string, must be 'pstrain' (plane strain). Other cases → error.
%   PROPMAT : 1×nmat struct array with per-material constants, e.g.:
%            .TYPECONST  ∈ {'ELAST','ELASTPLAS'}
%            .E          (Young’s modulus)
%            .nu         (Poisson’s ratio)
%            .dens       (density)
%            .ELASTPLAS.TYPE = 'LINEARISOTRO' (required for plastic vectorization)
%            .sigmay     (initial yield stress σ_y0)   [only for ELASTPLAS]
%            .Hmodul     (linear isotropic hardening H) [only for ELASTPLAS]
%   DATA : struct with runtime flags/aux; the routine sets defaults:
%            .StrainStressWith4Components = 1
%            .nstrain = 4
%            .ISPSTRAIN = 1
%            .MESH.ngaus must exist (used to expand elem → GP)
%
% OUTPUTS
%   MATPRO : struct with Gauss-point fields (length = nelem*ngausE):
%       .sigmay_0 (ngaus×1)  → initial yield stress per GP
%       .Hmodulus (ngaus×1)  → isotropic hardening modulus per GP
%       .mu       (ngaus×1)  → shear modulus μ = E / (2(1+ν))
%       .kappa    (ngaus×1)  → bulk modulus κ = E / (3(1−2ν))
%       .Celas    ( (ngaus*4)×4 ) → block-stacked 4×4 elastic matrices in Voigt
%       .identVOL ( (ngaus*4)×4 ) → repeated volumetric projector (I⊗I)
%       .identV   ( (ngaus*4)×4 ) → repeated deviatoric projector (with 0.5 on γ_xy)
%       .dens     (nelem×1)  → element density (copied from PROPMAT)
%   DATA : updated with fields described above
%
% METHOD (high level)
%   1) Check ndim and typePROBLEM → enforce plane strain (nstrain=4).
%   2) Normalize/extend MESH.MaterialType to match PROPMAT length.
%   3) For each material:
%        - Map its element set ELEMS to Gauss-point indices via small2large.
%        - Fill per-GP arrays: E, ν, (σ_y0, H) depending on TYPECONST.
%   4) Compute (μ, κ, λ) at each GP and assemble the 4×4 elastic tensor:
%        C = λ (I⊗I) + 2μ I_dev
%      implemented as:
%        identVOL = kron(ones(ngaus,1), [1 1 0 1]'*[1 1 0 1])   (block-repeat)
%        identV   = kron(ones(ngaus,1), diag([1, 1, 0.5, 1]))   (γ_xy scaling)
%        Celas    = λ*identVOL + 2μ*identV                      (vectorized)
%
% VOIGT ORDER & SCALING (plane strain, nstrain=4)
%   ε = [ε_xx, ε_yy, γ_xy, ε_zz]^T
%   Note: γ_xy is engineering shear (γ_xy = 2 ε_xy), hence identV(3,3)=0.5 so that
%         C_33 = 2μ (engineering shear modulus).
%
% UNITS
%   Use consistent units across E, σ_y0, H, ρ, etc. (e.g., MPa, m, kg/m^3).
%
% DEPENDENCIES
%   - DefaultField(DATA, fieldName, defaultVal) : set default flags
%   - small2large(elemSet, ngausE)               : map element indices to
%                                                 concatenated GP indices
%
% NOTES / PITFALLS
%   - For 'ELASTPLAS', only **parameters** (σ_y0, H) are set here; the plastic
%     return-mapping and algorithmic tangent are handled elsewhere.
%   - Plane stress and 3D not implemented here; extending requires adapting
%     nstrain, Voigt order, and the projection operators/closures for C.
%   - If MESH.MaterialType is empty, it defaults to ones(nelem,1).
%   - If NMAT in mesh > length(PROPMAT), last material is replicated.
%
% EXTENSIONS
%   - Add temperature-dependent properties by evaluating E(T), σ_y0(T), etc.
%   - Support anisotropy by providing per-material full 3D C and reducing
%     to plane strain/pstress consistently.
%   - Provide sparse block representation of C for memory efficiency.
% =========================================================================
% JAHO
% COMMENTS by CHATGPT5, 20-Oct-2025

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
MATPRO.dens = zeros(nelem,1) ;
%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
   % celas3D =PROPMAT(imat).ElasticityMatrix ; %
   % INVcelas3D = inv(celas3D) ;
    ELEMS = find(MESH.MaterialType == imat) ;    
%     switch typePROBLEM
%         case 'pstrain'
%             rowcol = [1 2 6] ;
%             celas = celas3D(rowcol,rowcol) ;
%         case 'pstress'
%             rowcol = [1 2 6] ;
%             celasINV3D = inv(celas3D) ;
%             celasINV = celasINV3D(rowcol,rowcol) ;
%             celas = inv(celasINV) ;
%         case '3D'
%             celas = celas3D ;
%     end
   % for eLOC=1:length(ELEMS)
    %    e = ELEMS(eLOC) ;
       % MATPRO.celasglo(:,:,e) = celas ;
        MATPRO.dens(ELEMS) = PROPMAT(imat).dens ;
       
    %end
end





DATA = DefaultField(DATA,'StrainStressWith4Components',1) ;  

DATA.nstrain = nstrain ; 
EXIST_DENS = 0 ;

DATA.ISPSTRAIN = 1; ;

%celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
%celasgloINV = zeros(6,6,nelem) ;  % Global array of compliance matrices (3D)
densGLO = ones(nelem,1) ;  % Global array of compliance matrices (3D)


NMAT_fe = length(unique(MESH.MaterialType));
nmat = length(PROPMAT) ;

if NMAT_fe > nmat
    for imat = nmat:NMAT_fe
        PROPMAT(imat) =  PROPMAT(nmat)   ;
    end
    nmat = NMAT_fe ;
end


if isempty(MESH.MaterialType)
    MESH.MaterialType = ones(nelem,1) ;
end


ngausE = DATA.MESH.ngaus ;
ngaus =  nelem*ngausE ; 
E_young = zeros(ngaus,1) ;  % Young's modulus
Hmodulus = zeros(ngaus,1) ; % Hardening parameters
PoissonCOEF = zeros(ngaus,1) ;  % Poisson's ratio
sigmay_0 = 1e20*ones(ngaus,1) ;  % Yield stress 
%dbstop('30')
for imat = 1:nmat
    % Load elasticity tensor
    MATLOC = PROPMAT(imat)  ;
    E = MATLOC.E ;
    nu = MATLOC.nu ;
    G = MATLOC.G ;
    
     
    elemMAT = find(MESH.MaterialType == imat) ;
    pointsMAT = small2large(elemMAT,ngausE) ; 
    
    
    switch PROPMAT(imat).TYPECONST
        case {'ELASTPLAS','ELAST'}
            E_young(pointsMAT) = PROPMAT(imat).E;
            PoissonCOEF(pointsMAT) = PROPMAT(imat).nu;
            switch  PROPMAT(imat).TYPECONST
                case 'ELASTPLAS'
                    switch   PROPMAT(imat).ELASTPLAS.TYPE
                        case 'LINEARISOTRO'
                            sigmay_0(pointsMAT) =  PROPMAT(imat).sigmay ;
                            Hmodulus(pointsMAT) = PROPMAT(imat).Hmodul;                            
                        otherwise
                            error('This problem is not amenable to vectorization (stress calculation)')
                    end
                    
            end
    end
    
    
     
    
    
    
    
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




% %%% Indices for obtaining Ctang
% [indI,indJ] = IndicesCtang(size(Celas,1),size(Celas,2)) ;
% MATPRO.indI = indI ;
% MATPRO.indJ = indJ ;
% 
% % Diagonal matrix (Elastic range)
%   MATPRO.CelasGLO =     ConvertCmatSparseMatrix(Celas,indI,indJ) ;
