function   [MESH,MATPRO,OPERHROM,DISP_CONDITIONS,DATAHROM,OTHER_data] ...
    =  PreProcessInputDataHOMOG_HROM(DATAHROM,DATAFE,MACRODEF,ECMdata,BasisU,BasisStwo,OTHER_INPUTS)
%--------------------------------------------------------------------------
% function [MESH,MATPRO,OPERHROM,DISP_CONDITIONS,DATAHROM,OTHER_data] ...
%   = PreProcessInputDataHOMOG_HROM(DATAHROM,DATAFE,MACRODEF,ECMdata,BasisU,BasisStwo,OTHER_INPUTS)
%--------------------------------------------------------------------------
% PURPOSE
%   Prepare all data structures required to run a large-strain FE-HROM with
%   periodic boundary conditions and PK1 homogenization. This routine
%   adapts the full-order FE inputs (materials, operators, mesh sizes) to
%   the reduced integration domain dictated by ECM (continuous or discrete),
%   checks BC consistency with the training setup, constructs the reduced
%   operators, and assembles homogenization weights.
%
%   It is a HROM-oriented modification of the FE preprocessor
%   `PreProcessInputDataDyn1.m`.
%
%--------------------------------------------------------------------------
% INPUT
%   DATAHROM     : Struct with HROM run parameters (e.g., time/path info).
%   DATAFE       : Struct produced in the FE offline stage containing:
%                  .FE_VARIABLES_NAMEstore → MAT file with MATPRO, MESH, OPERFE
%                  .MESH.ngaus_STRESS, .MESH.nstrain, .DOFr (training DOFs)
%   MACRODEF     : Definition of macroscopic kinematics (paths/steps), used
%                  by DispMACROtime to build MACROVAR (macro disp/grad fields).
%   ECMdata      : Empirical Cubature/Continuous ECM selection data:
%                  .setElements → selected elements (for continuous ECM)
%                  .setPoints   → selected Gauss points (for discrete ECM)
%                  .wRED        → reduced quadrature weights (stress points)
%                  .weightsHOMOG_PK1 → weights for PK1 homogenization
%   BasisU       : Reduced basis for displacement fluctuations (columns).
%   BasisStwo    : (Unused here) Basis for PK2 stresses; included for interface consistency.
%   OTHER_INPUTS : Optional struct for BC options (e.g., rigid-body flags). If
%                  omitted, defaults are assumed inside called routines.
%
%--------------------------------------------------------------------------
% OUTPUT
%   MESH              : Mesh struct (coordinates are replaced by relative/
%                       periodic coordinates returned by Periodic_BoundaryCOND_LARGE).
%   MATPRO            : Material database restricted to the reduced Gauss set:
%                       .celasglo → sliced to ECM points
%                       .dens     → sliced to ECM elements
%   OPERHROM          : Reduced operators for the HROM evaluation:
%                       .Ared        → map for fluctuation reconstruction (A*BasisU)
%                       .BstA        → reduced strain-disp operator at ECM points
%                       .wSTs        → reduced quadrature weights (stress)
%                       .IDENTITY_F  → identity F at ECM points
%                       .WEIGHTShomog→ weights for PK1 homogenization
%   DISP_CONDITIONS   : ROM boundary-condition struct:
%                       .MACROVAR → macro fields from DispMACROtime
%                       .DOFl     → 1:size(BasisU,2) (all ROM DOFs are free)
%                       .DOFr, .DOFm → empty (constraints enforced via periodicity)
%   DATAHROM          : Updated with mesh-related counts for the reduced grid:
%                       .MESH.nstrain, .MESH.ngaus_STRESS (copied from FE)
%                       .MESH.ngausT      → # of ECM stress points
%                       .MESH.ndofSTRESS  → ngausT * nstrain
%                       .MESH.ndof        → size(BasisU,2)
%   OTHER_data        : Reserved for auxiliary exports (currently empty).
%
%--------------------------------------------------------------------------
% METHOD / KEY STEPS
%   1) Load FE variables (MATPRO, MESH, OPERFE) from DATAFE.FE_VARIABLES_NAMEstore.
%   2) Restrict material arrays to ECM’s reduced integration set:
%        - Continuous ECM: map selected ELEMENTS to their Gauss rows using
%          `setPointsElement = setElements * ngaus_STRESS`, then expand to
%          stress components via `small2large(·, nstrain)`.
%        - Discrete ECM: select Gauss points directly via ECMdata.setPoints.
%        - Slice MATPRO.celasglo (constitutive tangents) with these indices.
%        - Slice MATPRO.dens by ECM elements.
%   3) Build periodic boundary conditions and relative coordinates:
%        [DISP_CONDITIONS_LOC, COORrel_LOC, ~] = Periodic_BoundaryCOND_LARGE(...)
%        Check that constrained DOFs match the training set (DATAFE.DOFr).
%        Replace MESH.COOR by COORrel_LOC' (reference for macro kinematics).
%   4) Macroscopic fields:
%        MACROVAR = DispMACROtime(MACRODEF, COORrel_LOC, DATAHROM);
%        Store in DISP_CONDITIONS.MACROVAR.
%   5) Define ROM DOF partition:
%        DOFl = 1:nmodesU (all ROM DOFs are free); DOFr, DOFm empty.
%   6) Reduced operators:
%        - BstRED = OPERFE.BstA * BasisU.
%        - If continuous ECM: interpolate BstRED to ECM Gauss points via
%          InterpolationGaussVariablesECM(BstRED, ECMdata, ngaus_STRESS, nF).
%          Else: pick rows by `small2large(ECMdata.setPoints, nF)`.
%        - Set OPERHROM.Ared = DISP_CONDITIONS_LOC.A * BasisU (fluctuation reconstruction).
%        - Set weight vectors and identity deformation entries at ECM points.
%   7) Mesh counters (reduced grid): ngausT, ndofSTRESS, ndof.
%   8) Homogenization weights for PK1 stresses:
%        OPERHROM.WEIGHTShomog = ECMdata.weightsHOMOG_PK1.
%
%--------------------------------------------------------------------------
% ASSUMPTIONS / REQUIREMENTS
%   * Periodic BCs were used in training; the current constrained DOFs must
%     match DATAFE.DOFr. Otherwise, execution stops with an error.
%   * For continuous ECM, material properties are assumed constant across
%     Gauss points of a given element for slicing consistency.
%   * BasisU must be column-orthonormal (or well-conditioned) to avoid loss
%     of accuracy in BstRED and Ared.
%
%--------------------------------------------------------------------------
% NOTES (Implementation)
%   * `small2large(idxs, ncomp)` expands point indices to block-row indices
%     accounting for ncomp stress components per Gauss point.
%   * `InterpolationGaussVariablesECM` maps quantities stored per-element/
%     per-Gauss in the full grid to the ECM point set (continuous ECM).
%   * `OPERFE.IDENTITY_F` is reshaped/sliced to the ECM point set so that
%     F = I + (·) can be formed consistently at reduced quadrature locations.
%   * `BasisStwo` is not used here but retained for interface compatibility
%     with pipelines that may require it downstream.
%
%--------------------------------------------------------------------------
% AUTHOR
%   Joaquín A. Hernández Ortega (UPC/CIMNE) — jhortega@cimne.upc.edu
%
% DATE
%   Comments by ChatGPT-6-October-2025, Barcelona
%
%--------------------------------------------------------------------------
% GLOSARIO (rápido, ES)
%   * continuous ECM: ECM continuo (selección por elementos + interpolación)
%   * discrete ECM: ECM discreto (selección directa de puntos de Gauss)
%   * slicing: restricción/selección de filas (recorte de arreglos)
%   * homogenization weights: pesos de homogenización (PK1)
%--------------------------------------------------------------------------


%
OTHER_data =[] ;
if nargin == 0
    load('tmp1.mat')
elseif nargin == 8
    OTHER_INPUTS = []  ; 
end


% ---------------
% RECOVER MATPRO
%---------------
load(DATAFE.FE_VARIABLES_NAMEstore,'MATPRO','OTHER_output','MESH','OPERFE') ;
if  isempty(ECMdata.setPoints)
    % The continuous ECM has been used 
    % We assume that, for a given element, the elastic properties are the
    % same for all the Gauss points of such elements 
    setPointsElement = ECMdata.setElements*DATAFE.MESH.ngaus_STRESS; % Last Gauss point  of an element
    % For instance, if ngaus_STRESS = 9, and setElement = 1, then 
    % setPointsElement=9 
    setIndicesLOC = small2large(setPointsElement,DATAFE.MESH.nstrain) ;
    MATPRO.celasglo = MATPRO.celasglo(setIndicesLOC,:) ;
    
else
    setIndices = small2large(ECMdata.setPoints,DATAFE.MESH.nstrain) ;
    MATPRO.celasglo = MATPRO.celasglo(setIndices,:) ;
    
end
MATPRO.dens = MATPRO.dens(ECMdata.setElements,:) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RB_MOTION',[]) ; 
% 
% [DOFrFE,dRfe] = DirichletCONDtime(DIRICHLET,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ;
% III =  setdiff(DOFrFE,DATAFE.DOFr) ;
% if ~isempty(III)
%     error('The constrained DOFs should be the same that the one used for training')
% end


[DISP_CONDITIONS_LOC,COORrel_LOC,DATA_LOC] = Periodic_BoundaryCOND_LARGE(MESH.COOR,MESH.CN,MESH.CNb,DATAFE) ;
 III =  setdiff(DISP_CONDITIONS_LOC.DOFr,DATAFE.DOFr) ;
 if ~isempty(III)
     error('The constrained DOFs should be the same that the one used for training')
 end

MESH.COOR = COORrel_LOC' ; 
DATAHROM.ngausTotalLOC = length(ECMdata.wRED) ;
MACROVAR = DispMACROtime(MACRODEF,COORrel_LOC,DATAHROM) ; 
DISP_CONDITIONS.MACROVAR = MACROVAR ; 
DISP_CONDITIONS.DOFr = [] ; 
DISP_CONDITIONS.DOFl = 1:size(BasisU,2) ; 
DISP_CONDITIONS.DOFm = [] ; 


% OTHER_data.DOFlFE = DISP_CONDITIONS_LOC.DOFl ; 
% OTHER_data.DOFrFE = DISP_CONDITIONS_LOC.DOFr ; 
 
%%%%%%
% OPERFE.BstA = OPERFE.BstA(:,:)*BasisU  ; 





%%%%%%%% DEFINITION OF BasisUall
% --------------------------------------
%BasisU_r = dRfe.U;   % Constrained part
%nmodesUr = size(BasisU_r,2) ;
%nmodesUall = nmodesUr + nmodesUl ;
%BasisUall = zeros(DATAFE.MESH.ndof,nmodesUall) ;
%COLS  = 1:nmodesUl ;
%BasisUall(DOFlFE,COLS) = BasisU ;
%COLS  = nmodesUl+1:nmodesUall ;
%BasisUall(DOFrFE,COLS) = BasisU_r ;
%DATAHROM.MESH.ndof = nmodesUall;
% ------------------------------------------
% DIRICHLET CONDITIONS REDUCED-ORDER MODEL
% -------------------------------------------------
%nmodesUl = size(BasisU,2) ;
%DISP_CONDITIONS.DOFl = 1:nmodesUl ;
%DISP_CONDITIONS.DOFr = nmodesUl+1:nmodesUall;
%dR.U = eye(nmodesUr);
%dR.a = dRfe.a ;
%DISP_CONDITIONS.dR = dR  ;
%%% OPERATORS
% -----------
BstRED = OPERFE.BstA*BasisU ;
OPERHROM.Ared = DISP_CONDITIONS_LOC.A*BasisU ; % For reconstruction of fluctuations 
nF = DATAFE.MESH.ndim^2 ;
if isempty(ECMdata.setPoints)
    % Interpolation  (  we are using the Continuous ECM )
    % --------------------------------------------------------------------  
    OPERHROM.BstA =  InterpolationGaussVariablesECM(BstRED,ECMdata,DATAFE.MESH.ngaus_STRESS,nF) ;
    setIndices = small2large(setPointsElement,nF) ;
else
    setIndices = small2large(ECMdata.setPoints,nF) ;
    OPERHROM.BstA = BstRED(setIndices,:) ;
end
OPERHROM.wSTs =  ECMdata.wRED;
OPERHROM.IDENTITY_F =  OPERFE.IDENTITY_F(setIndices);
DATAHROM.MESH.nstrain = DATAFE.MESH.nstrain ;
DATAHROM.MESH.ngaus_STRESS = DATAFE.MESH.ngaus_STRESS;

DATAHROM.MESH.ngausT =  length(OPERHROM.wSTs) ; % Total number of Gauss points
DATAHROM.MESH.ndofSTRESS =   DATAHROM.MESH.ngausT*DATAHROM.MESH.nstrain ; % Total number of Gauss points
DATAHROM.MESH.ndof = size(BasisU,2);


OPERHROM.WEIGHTShomog = ECMdata.weightsHOMOG_PK1 ; 
 % HOMOGENIZATION PK1 STRESSES
 
