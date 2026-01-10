function [BASIS_DEF,MESH,DATA,DATAoffline,Kstiff,MdomCHOL,MintfINV_chol] ...
    = DeformModesNONL_INVAR(OPERFE,MESH,DATA,...
    SNAPdisp,DATAcommon,DATAoffline,Vintf,PhiRB,Mdom,Mintf,MATPRO,SNAPdisp_AUX)
%--------------------------------------------------------------------------
% FUNCTION: DeformModesNONL_INVAR
%
% PURPOSE:
%   Constructs a set of deformational basis modes for a given subdomain 
%   within a nonlinear multiscale Reduced Order Model (ROM), such as EIFEM.
%   The approach distinguishes between:
%     - Elastic (invariant) deformation patterns extracted from training tests,
%     - Inelastic (complementary) deformation modes from more general tests.
%   These are separated and orthogonalized with respect to rigid body motions.
%
%   The method relies on SVD-based extraction of dominant strain modes,
%   uses Cholesky-based mass inner products, and separates contributions
%   from multiple domains (if auxiliary snapshots are provided).
%
% INPUT ARGUMENTS:
%   - OPERFE      : Structure containing FE operators (B-matrix, Gauss weights, etc.)
%   - MESH        : Mesh structure for the current subdomain
%   - DATA        : Configuration and solver options
%   - SNAPdisp    : Cell array of displacement snapshots (ordered: elastic first, then inelastic)
%   - DATAcommon  : Shared/common data across subdomains (e.g., tolerances, projections)
%   - DATAoffline : Configuration for offline phase (e.g., snapshot scaling, index tracking)
%   - Vintf       : Interface basis functions (e.g., projections on rigid-body DOFs)
%   - PhiRB       : Matrix of rigid body displacement modes
%   - Mdom        : Domain mass matrix (for inner products and orthogonalizations)
%   - Mintf       : Interface mass matrix (used in WSVDT and optional orthonormalization)
%   - MATPRO      : Structure with material properties (`celasglo` or `Celas`)
%   - SNAPdisp_AUX: Optional cell array of additional displacement snapshots 
%                   from auxiliary domains (used for multiblock training)
%
% OUTPUT ARGUMENTS:
%   - BASIS_DEF       : Structure containing:
%                        * UpsilonDEF     : Full set of orthonormal strain modes
%                        * S_upsilonDEF   : Singular values associated with strain modes
%                        * V_upsilonDEF   : Right singular vectors from SVD
%                        * PhiDEFbs       : Elastic (basic) deformation modes
%                        * PhiDEFcomp     : Inelastic (complementary) modes (if any)
%   - MESH            : Possibly updated mesh structure (e.g., with face DOFs)
%   - DATA            : Updated configuration, including mode-related fields
%   - DATAoffline     : Updated offline configuration (e.g., index bookkeeping)
%   - Kstiff          : Stiffness matrix of the domain (either assembled or recovered)
%   - MdomCHOL        : Cholesky factor of the domain mass matrix
%   - MintfINV_chol   : Inverse Cholesky factor of the interface mass matrix
%
% ALGORITHM OVERVIEW:
%   STEP 1: Split the displacement snapshots into elastic and inelastic parts
%   STEP 2: Rescale basic (elastic) and complementary (inelastic) snapshots
%           possibly incorporating auxiliary domains
%   STEP 3: Recover stiffness matrix based on material tensor and FE operators
%   STEP 4: Compute elastic deformation modes using constrained SVD
%   STEP 5: Build complementary modes by projecting out basic modes,
%           followed by WSVDT-based orthogonalization
%   STEP 6: Collect strain modes via SVD on inelastic snapshots,
%           truncate to the specified number of dominant modes
%
% FEATURES:
%   - Compatible with both single-domain and multi-domain training setups
%   - Elastic/inelastic mode separation
%   - Orthogonality with respect to both mass matrix and rigid body modes
%   - Fully SVD-driven construction
%
%
% RELATED FUNCTIONS:
%   - InvariantModes_basicNC
%   - ComplementaryModes_NOTESTS
%   - StrainModesNonlinearSVD_EIFEM
%   - WSVDT
%   - SprojDEF_operator
%   - PlotModesDEF_SubdomainLevel
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC / CIMNE)
%   First version: 03-Apr-2024, Barcelona
%   Comments and documentation update: ChatGPT-4, 1-JUN-2025
%--------------------------------------------------------------------------
 
%% JAHO, 3-Apr-2024
%  Barcelona, Sandwichez, C/Mallorca
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% Determination of deformational modes for the case in which we are
% interesed in computing "invariant" modes, such as "beam" modes
% See training stage for elastic range
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE/Train_Q8_48tests.m

if nargin == 0
    load('tmp1.mat')
end

% STEP 1
% --------------------------------------------------------------------------
% Nodes and DOFs of the interface boundaries
%[faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESH] =   GeometricVarDOMAINS(OPERFE,MESH,DATA,DATAcommon);

% STEP 2
%%%%%%%%---------------------------
% Basic tests/ Complementary (inelastic tests) tests
% --------------------------------------
ntestsCOMPL = length(DATAoffline.AdditionalTests) ;  % THESE ARE THE INELASTIC TESTS
ntestsALL = length(SNAPdisp) ;
INDEX_BASIC_TESTS = 1:(ntestsALL-ntestsCOMPL) ;   % THESE ARE THE INDICES OF THE ELASTIC MDES
INDEX_COMPL_TESTS = (INDEX_BASIC_TESTS(end)+1):ntestsALL;
DATAoffline.INDEX_BASIC_TESTS = INDEX_BASIC_TESTS;
DATAoffline.INDEX_COMPL_TESTS = INDEX_COMPL_TESTS;

%  Scales displacement snapshot data (e.g., from training simulations) 
%   according to predefined influence weights for each domain, allowing 
%   the user to modulate the contribution of different domains in the 
%   computation of deformation modes (e.g., in EIFEM or ROM settings).
[SNAPbasic,SNAPcompl,DATAoffline] = ScalingSnapshotsDISP(SNAPdisp,SNAPdisp_AUX,INDEX_BASIC_TESTS,DATAoffline,INDEX_COMPL_TESTS);


% Stiffness matrix
if isfield(MATPRO,'Celas')
    MATPRO.celasglo = MATPRO.Celas ;
end

if size(MATPRO.celasglo,1) == size(OPERFE.Bst,1)
    Kstiff = StiffMatrixRecover(MATPRO,OPERFE)  ;
else
    Kstiff = StiffMatrixRecoverLARGE(MATPRO,OPERFE,MESH,DATA)  ;
    
end


% step 3) INVARIANT MODES (LINEAR)
% -------------------------
[PhiDEFbs_INV,MdomCHOL,DATAoffline] = InvariantModes_basicNC(SNAPbasic,DATAcommon,PhiRB,Mdom,DATAoffline,MESH,DATA,Kstiff,Mintf) ;

% step 4) COMPLEMENTARY BASIC MODES (STILL LINEAR) --> bASIC MODES


[PhiDEFbs,MintfINV_chol,DATA] = ComplementaryModes_NOTESTS...
    (PhiDEFbs_INV,DATAcommon,SNAPbasic,MdomCHOL,MESH,PhiRB,Vintf,Mintf,DATA,MATPRO,OPERFE,Mdom,Kstiff,DATAoffline) ;




%%% STEP 5. PURGING RIGID BODY DISPLACEMENTS, complementary tests. Basis
%%% PhiDEFcomp

disp('-----------------------------------------------')
disp('Determining nonlinear strain modes ')
disp('-------------------------------------------------')


[UpsilonDEF,S_upsilonDEF,V_upsilonDEF] = StrainModesNonlinearSVD_EIFEM(SNAPcompl,DATAoffline,Mdom,PhiDEFbs,PhiRB,MdomCHOL,MESH)  ;



DATAoffline = DefaultField(DATAoffline,'NUMBER_OF_DEFORMATIONAL_MODES',[]) ;

if  ~isempty(DATAoffline.NUMBER_OF_DEFORMATIONAL_MODES)
    nmodes = length(S_upsilonDEF) ;
    nmodes = min(nmodes,DATAoffline.NUMBER_OF_DEFORMATIONAL_MODES) ;
else
    nmodes = length(S_upsilonDEF) ;
end


BASIS_DEF.UpsilonDEF = UpsilonDEF(:,1:nmodes) ;
BASIS_DEF.S_upsilonDEF = S_upsilonDEF(1:nmodes) ;
BASIS_DEF.V_upsilonDEF = V_upsilonDEF(:,1:nmodes) ;



DATA.NAME_BASE = 'AllDEFORMATIONAL';
PlotModesDEF_SubdomainLevel(DATA,UpsilonDEF,MESH);

% %% OUT OF CURIOSITY, HOW DOES THE SVD OF THE BOUNDARY ENTRIES LOOK LIKE ?
%                 b = MESH.faceDOFSall;
%
% Ub  = bsxfun(@times,UpsilonDEF(b,:)',S_upsilonDEF)' ;
% DATALOCaaa.TOL = 1e-4 ;
% [ Ub_b,S_b,V_b] = WSVDT( Ub,Mintf,DATALOCaaa) ;
%


% ORTHOGONAL COMPLEMENT


ncomp  = size(BASIS_DEF.UpsilonDEF,2)-size(PhiDEFbs,2) ;

% if ncomp >1 && ~isempty(SNAPcompl)

if ncomp >0&& ~isempty(SNAPcompl)  % Amended on Sept-18th-2024
    
    PhiDEFcomp = SprojDEF_operator(PhiDEFbs,Mdom, BASIS_DEF.UpsilonDEF) ;
    DATALOC.TOL = 1e-6 ;
    DATALOC.Mchol = MdomCHOL ;
    [ PhiDEFcomp,S,~,~] = WSVDT( PhiDEFcomp,[],DATALOC) ;
    
    PhiDEFcomp = PhiDEFcomp(:,1:ncomp)  ;
    
    DATA.NAME_BASE = 'Complementary_orth';
    PlotModesDEF_SubdomainLevel(DATA,PhiDEFcomp,MESH);
    
    
else
    disp(['There are no complementary (inelastic modes ) here, only elastic'])
    PhiDEFcomp = [] ;
end



BASIS_DEF.PhiDEFbs = PhiDEFbs ;
BASIS_DEF.PhiDEFcomp = PhiDEFcomp ;




