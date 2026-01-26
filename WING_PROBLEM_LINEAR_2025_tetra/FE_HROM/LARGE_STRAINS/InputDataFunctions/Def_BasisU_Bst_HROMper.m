function [DISP_CONDITIONS_HROM,OPERHROM,OTHER_data,BasisUall] = ...
    Def_BasisU_Bst_HROMper(BasisU,DATAFE,OTHER_INPUTS,DATAHROM,OTHER_data,OPERFE,ECMdata,DISP_CONDITIONS)
%==========================================================================
% Def_BasisU_Bst_HROMper  —  (affine Dirichlet BCs)
%
% PURPOSE
%   Build the **global ROM basis** and **hyper-reduced stress operator** when
%   essential BCs are imposed through an **affine map**. In this setting, the
%   full displacement lives on a manifold and is represented with a basis that
%   already satisfies the Dirichlet constraints.
%
%   The *global* reduced basis is (see screenshot / Eq. (257) in the slides):
%
%         Φ̄ = [A Φ   0]  +  [ 0  0
%                              0  D^U ]
%
%   where:
%     • Φ   : free-DOF basis (columns on DOFl),
%     • A   : affine inserter mapping free DOFs into the full vector,
%     • D^U : boundary/Dirichlet columns spanning prescribed DOFs (on DOFr).
%
%   “Standard” BCs are recovered by choosing A(l,:) = I and A(r,:) = 0, i.e.,
%   Φ̄ = [Φ  0] + [0 0; 0 D^U]. This matches the weak-form view where trial
%   functions satisfy essential BCs and unknowns are the free block, with the
%   usual block elimination Kll dl = Fl − Klr dr. 
%
% WHAT THIS ROUTINE DOES
%   1) Forms BasisUall by placing A*BasisU in *all* rows (affine injection) and
%      appending the boundary columns D^U on DOFr, exactly as Φ̄ above.
%   2) Defines the ROM Dirichlet data in reduced coordinates:
%         y = [ y_l ; y_r ],  with d_R(t) = D^U * a(t),
%      i.e., dR.U = I (on the y_r block) and dR.a = DISP_CONDITIONS.dR.a.
%   3) Builds the hyper-reduced operator on the ECM set:
%         BstRED = Bst * BasisUall,
%         OPERHROM.Bst = BstRED(restricted to ECM points).
%      This follows the slides’ vectorized assembly F̂ = Bᵀ ω S and its ECM
%      restriction to selected Gauss points. 
%   4) Carries over IDENTITY_F at ECM points when the formulation uses F.
%
% KEY DETAILS / CONSISTENCY NOTES
%   • Partition in ROM space:
%         DISP_CONDITIONS.DOFl = 1 : nRED,
%         DISP_CONDITIONS.DOFr = nRED+1 : ndof_ROM,
%     with ndof_ROM = nRED + size(D^U,2).
%   • Choice of strain stack length nF:
%       – small-strain w/o F: nF = nstrain (Voigt),
%       – otherwise (finite strain or F-based small strain): nF = ndim^2.
%   • ECM mapping:
%       – If ECMdata.setPoints is empty → continuous ECM interpolation.
%       – Else → discrete selection via small2large(…, nF) on BstRED rows.
%   • Theory tie-ins:
%       – Essential BCs enforced in the trial space (trial functions satisfy u=g),
%         while natural BCs arise from the weak form. :contentReference[oaicite:2]{index=2}
%       – Block decomposition and elimination for dl align with the classic
%         Kll/Klr form. :contentReference[oaicite:3]{index=3}
%       – Vectorized operators (L,B,N,ω) justify precomputing Bst and then
%         restricting to ECM points. :contentReference[oaicite:4]{index=4}
%
% INPUTS
%   BasisU          : Free-DOF basis Φ (on DOFl in training coordinates).
%   DATAFE, OTHER_INPUTS, DATAHROM, OTHER_data, OPERFE, ECMdata : usual FE/HROM
%                        structures; ECMdata may carry setPoints / setElements.
%   DISP_CONDITIONS : Must contain fields A (affine inserter), DOFr, dR.U, dR.a.
%
% OUTPUTS
%   DISP_CONDITIONS_HROM : Dirichlet info in ROM coordinates (dR.U = I on y_r).
%   OPERHROM.Bst         : Stress operator on ECM points in ROM coordinates.
%   OPERHROM.IDENTITY_F  : Identity contribution at ECM points (if available).
%   OTHER_data.BasisU_r  : Stored Dirichlet columns (D^U).
%   BasisUall            : Global ROM basis Φ̄ consistent with affine BCs.
%
% PITFALLS
%   • A, DOFr, and D^U must be mutually consistent and match the training split.
%   • nF must match the chosen kinematics branch; otherwise small2large indices
%     will be misaligned.
%==========================================================================
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
% JAHO, 8-Oct-2025, 19:23, Starbucks, Paseo de Gracia, Barcelona
% --------------------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

%%%%%%%% DEFINITION OF BasisUall
% --------------------------------------
BasisU_r = full(DISP_CONDITIONS.dR.U);   % Constrained part
nmodesUr = size(BasisU_r,2) ;
nmodesUl = size(BasisU,2) ;
nmodesUall = nmodesUr + nmodesUl ;
BasisUall = zeros(DATAFE.MESH.ndof,nmodesUall) ;
COLS  = 1:nmodesUl ;
BasisUall(:,COLS) = DISP_CONDITIONS.A*BasisU ;
%DATA_BasisUAll.DOFlFE = DOFlFE ;
%DATA_BasisUAll.COLS_DOFlFE = DOFlFE ;
%------------
COLS  = nmodesUl+1:nmodesUall ;
BasisUall(DISP_CONDITIONS.DOFr,COLS) = BasisU_r ;
DATAHROM.MESH.ndof = DATAHROM.nREDcoord+size(BasisU_r,2);
%------------
%DATA_BasisUAll.DOFrFE = DOFrFE ;
%DATA_BasisUAll.COLS_DOFrFE = COLS ;
%DATA_BasisUAll.BasisU_r = BasisU_r;
OTHER_data.BasisU_r = BasisU_r;
% ------------------------------------------
% DIRICHLET CONDITIONS REDUCED-ORDER MODEL
% -------------------------------------------------OTHER_data
DISP_CONDITIONS_HROM.DOFl = 1:DATAHROM.nREDcoord ;
DISP_CONDITIONS_HROM.DOFr = (DATAHROM.nREDcoord+1):DATAHROM.MESH.ndof;
dR.U = eye(nmodesUr);
dR.a = DISP_CONDITIONS.dR.a ;
DISP_CONDITIONS_HROM.dR = dR  ;
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