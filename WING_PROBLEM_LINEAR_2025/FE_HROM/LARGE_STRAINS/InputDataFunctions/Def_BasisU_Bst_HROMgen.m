function [DISP_CONDITIONS,OPERHROM,OTHER_data,BasisUall] = ...
    Def_BasisU_Bst_HROMgen(DOFrFE,dRfe,BasisU,DATAFE,OTHER_INPUTS,DATAHROM,OTHER_data,OPERFE,ECMdata,DOFlFE,DISP_CONDITIONS)
%==========================================================================
% Def_BasisU_Bst_HROMgen
%
% PURPOSE
%   Build the **global reduced basis** and the **hyper-reduced stress operator**
%   for a nonlinear HROM on a τ(q) manifold, consistent with **affine Dirichlet
%   (master–slave) boundary conditions**:
%       • Assemble BasisUall by stacking:
%           – columns for FREE DOFs (BasisU on DOFlFE),
%           – columns for PRESCRIBED/SLAVE DOFs taken from Dirichlet modes (dRfe.U).
%       • Map the full FE stress operator Bst to:
%           – the reduced coordinates via BstRED = Bst * BasisUall,
%           – the **selected Gauss points** (ECM set) to form OPERHROM.Bst.
%       • Provide DISP_CONDITIONS for the ROM partition (free/prescribed in
%         reduced coords) and pack auxiliary data in OTHER_data.
%
% AFFINE DIRICHLET, PARTITION, AND ROM SIZES
%   • Full FE space is split into DOFlFE (free) and DOFrFE (prescribed).
%   • The global reduced vector is partitioned as:
%         y = [ y_l ; y_r ],   with
%         y_l ∈ R^{nRED}  (free manifold coordinates),
%         y_r ∈ R^{nmodesUr} (amplitudes spanning constrained DOFs).
%   • BasisUall collects both blocks so that the full displacement is
%         d ≈ BasisUall * y,
%     with columns restricted to their respective index sets (DOFlFE/DOFrFE).
%   • DISP_CONDITIONS.DOFl = 1:nRED,  DISP_CONDITIONS.DOFr = nRED+1:ndof_ROM,
%     and Dirichlet time law is carried by dR (identity on y_r with amplitudes a(t)).
%
% CHOICE OF STRAIN VECTOR LENGTH nF
%   • For **small strains** without using F (deformation gradient) explicitly,
%     set nF = DATAFE.MESH.nstrain (Voigt size).
%   • Otherwise (finite strains or small strains using F), use nF = ndim^2,
%     matching the vectorized F-based formulation used downstream.
%
% ECM MAPPING (GAUSS-POINT SELECTION)
%   • If ECMdata.setPoints is **empty** → continuous ECM path:
%       OPERHROM.Bst = InterpolationGaussVariablesECM(BstRED, …)
%       OPERHROM.IDENTITY_F is extracted consistently for those points.
%   • If ECMdata.setPoints is **given** → discrete ECM path:
%       Restrict rows of BstRED (and IDENTITY_F if present) to the indices
%       corresponding to the selected Gauss points, expanded by nF via
%       small2large(…, nF).
%   • OPERHROM.wSTs (weights) are set elsewhere; here we only build operators.
%
% SPECIAL CASE: MISMATCHED TRAINING/TEST DOFr
%   • If OTHER_INPUTS.BasisU_allDOFS is provided (e.g., different boundary
%     partitions between training and runtime), re-estimate the free-part
%     basis via an SVD on the submatrix restricted to DOFlFE before building
%     BasisUall.
%
% OUTPUTS
%   BasisUall        : Global reduced basis (full-space ndof × (nRED + nmodesUr)).
%   DISP_CONDITIONS  : ROM partition and Dirichlet amplitudes (dR.U = I, dR.a = dRfe.a).
%   OPERHROM.Bst     : Hyper-reduced stress operator on ECM points.
%   OPERHROM.IDENTITY_F : Identity contribution at ECM points (if available).
%   OTHER_data.BasisU_r : Stored constrained/boundary columns used in BasisUall.
%
% PRACTICAL NOTES / PITFALLS
%   • Ensure DOFrFE matches the offline/training partition; otherwise the
%     Dirichlet block of BasisUall (and thus dR) becomes inconsistent.
%   • The mapping small2large(…, nF) expands Gauss-point indices to stacked
%     strain components; nF must match the chosen kinematics branch.
%   • If OPERFE.IDENTITY_F is empty (e.g., pure small-strain with no F),
%     leave OPERHROM.IDENTITY_F empty as well.
%==========================================================================
% JAHO, 8-Oct-2025, 19:19, Starbucks, Paseo de Gracia, Barcelona
% Comments by ChatGPT5

%%%%%%%% DEFINITION OF BasisUall
% --------------------------------------
BasisU_r = full(dRfe.U);   % Constrained part
nmodesUr = size(BasisU_r,2) ;
nmodesUl = size(BasisU,2) ;
nmodesUall = nmodesUr + nmodesUl ;
BasisUall = zeros(DATAFE.MESH.ndof,nmodesUall) ;
COLS  = 1:nmodesUl ;

OTHER_INPUTS = DefaultField(OTHER_INPUTS,'BasisU_allDOFS',[]); 

if ~isempty(OTHER_INPUTS.BasisU_allDOFS)
    % Situations in which the DOFr for training and testing are not the
    % same
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/02_HROMstand.mlx
    BasisU = SVDT(OTHER_INPUTS.BasisU_allDOFS(DOFlFE,:)) ; 
end



BasisUall(DOFlFE,COLS) = BasisU ;
%DATA_BasisUAll.DOFlFE = DOFlFE ;
%DATA_BasisUAll.COLS_DOFlFE = DOFlFE ;
%------------
COLS  = nmodesUl+1:nmodesUall ;
BasisUall(DOFrFE,COLS) = BasisU_r ;
DATAHROM.MESH.ndof = DATAHROM.nREDcoord+size(BasisU_r,2);
%------------
%DATA_BasisUAll.DOFrFE = DOFrFE ;
%DATA_BasisUAll.COLS_DOFrFE = COLS ;
%DATA_BasisUAll.BasisU_r = BasisU_r;
OTHER_data.BasisU_r = BasisU_r;
% ------------------------------------------
% DIRICHLET CONDITIONS REDUCED-ORDER MODEL
% -------------------------------------------------OTHER_data
DISP_CONDITIONS.DOFl = 1:DATAHROM.nREDcoord ;
DISP_CONDITIONS.DOFr = (DATAHROM.nREDcoord+1):DATAHROM.MESH.ndof;
dR.U = eye(nmodesUr);
dR.a = dRfe.a ;
DISP_CONDITIONS.dR = dR  ;
%%% OPERATORS
% -----------
BstRED = OPERFE.Bst*BasisUall ;
DATAFE = DefaultField(DATAFE,'SMALL_STRAIN_KINEMATICS',0) ; %  =

DATAFE = DefaultField(DATAFE,'NO_USE_Deformation_gradient_in_Small_Strains',1) ; %  =

if DATAFE.SMALL_STRAIN_KINEMATICS == 1 && DATAFE.NO_USE_Deformation_gradient_in_Small_Strains == 1
    nF = DATAFE.MESH.nstrain ;
else
    nF = DATAFE.MESH.ndim^2 ;
end


if isempty(ECMdata.setPoints)
    % Interpolation  (  we are using the Continuous ECM )
    % --------------------------------------------------------------------
    OPERHROM.Bst =  InterpolationGaussVariablesECM(BstRED,ECMdata,DATAFE.MESH.ngaus_STRESS,nF) ;
    setIndices = small2large(1:length(ECMdata.setElements),nF) ;  % This operation should be written more efficiently !!! 3-DEc-2021
    OPERHROM.IDENTITY_F =  OPERFE.IDENTITY_F(setIndices);
else
    setIndices = small2large(ECMdata.setPoints,nF) ;
    OPERHROM.Bst = BstRED(setIndices,:) ;
    if  ~isempty(OPERFE.IDENTITY_F)
        OPERHROM.IDENTITY_F =  OPERFE.IDENTITY_F(setIndices); % This operation should be removed !!! 3-DEc-2021
    else
        OPERHROM.IDENTITY_F = [] ;
    end
end