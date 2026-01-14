function [PhiLIN,PhiNON,qMASTER,qNON,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output,DATAoffline]...
    = Determine_qinf_qSUP_1NODE(CASES,NAMEsnap_base,DATAoffline,NAME_BASE,DATA_interp)
%==========================================================================
% Determine_qinf_qSUP_1NODE is a modification of
% Determine_qinf_qsup_LocMaxDisp, described below
% Its goal is to use as latent coordinate the displacement of a single node
% JAHO, 17-Nov-2025, Balmes 185, Barcelona

% Determine_qinf_qsup_LocMaxDisp is a modification of
% Determine_qinf_qsup_1SVDd, described in detail below. 
% The goal of the modification is to adapt the function so that the latent
% variable is the location of the maximum of the displacement in a given
% direction, and in a given surface
% JAHO, 16-Nov-2025, Sunday, 10:52. Balmes 185, BArcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/19_SHIFT_COMPRESS.mlx

%
% File: Determine_qinf_qsup_1SVDd.m
%
% Purpose
% -------
%   Constructs the reduced-order basis and generalized coordinates (qINF, qSUP)
%   from multiple sets of displacement snapshots, properly handling
%   non-homogeneous Dirichlet boundary conditions.  
%
%   This function is a refined version of `Determine_qinf_qsup_1SVD`, adapted
%   to problems where some prescribed displacements are **non-zero**. It ensures
%   the consistent extraction of linear (Φₗᵢₙ) and nonlinear (Φₙₒₙ) modes,
%   suitable for manifold-based Reduced Order Models (ROMs).
%
% Context
% -------
%   Used in nonlinear manifold-ROM workflows for both small- and
%   large-strain problems, where the reduced kinematic representation is
%   defined by:
%
%       u(x, t) ≈ Φₗᵢₙ·qₗᵢₙ(t) + Φₙₒₙ·qₙₒₙ(t)
%
%   The extracted qINF (scalar) and qSUP (vector) generalized coordinates
%   characterize the parametric evolution of the nonlinear manifold.
%
% Function Signature
% ------------------
%   [PhiLIN, PhiNON, qINF, qSUP, UU, SS, VV, ...
%    MESH, DATA, DOFl, DOFr, OPERFE, MATPRO, DISP_CONDITIONS, OTHER_output] = ...
%       Determine_qinf_qsup_1SVDd(CASES, NAMEsnap_base, DATAoffline)
%
% Inputs
% ------
%   CASES         : Array of integer IDs of training folders to process.
%   NAMEsnap_base : Base path to the snapshot directories.
%   DATAoffline   : Struct with offline parameters and tolerances:
%                     • errorDISP : SVD truncation tolerance for displacements.
%
% Outputs
% -------
%   PhiLIN          : [nDOF×1] Linear mode (first left singular vector).
%   PhiNON          : [nDOF×r] Nonlinear basis (remaining left singular vectors).
%   qINF            : [1×n_snap] Reduced coordinate on ΦLIN.
%   qSUP            : [r×n_snap] Reduced coordinates on ΦNON.
%   UU, SS, VV      : SVD of qSUP (used for further reduction or fitting).
%   MESH, DATA      : FE mesh and snapshot metadata (from first case).
%   DOFl, DOFr      : Indices of free and constrained DOFs.
%   OPERFE          : Finite element operators.
%   MATPRO          : Material properties.
%   DISP_CONDITIONS : Dirichlet boundary condition structure.
%   OTHER_output    : Auxiliary info (Dirichlet data, plotting basis, etc.).
%
% Method
% -------
%   1. Load displacement snapshots for all CASES.
%   2. Concatenate snapshots and reconstruct U*S*Vᵀ blocks.
%   3. Compute an overall SVD (for both constrained and unconstrained DOFs)
%      to obtain a visualization basis (OTHER_output.Phi_To_Plot).
%   4. Restrict each snapshot to the unconstrained DOFs (DOFl).
%   5. Apply weighted SVD (SRSVD) to extract Φ and singular values.
%   6. Split Φ = [Φₗᵢₙ Φₙₒₙ].
%   7. Compute reduced coordinates:
%         qINF = Φₗᵢₙᵀ·U,   qSUP = Φₙₒₙᵀ·U
%   8. Remove duplicates in qINF (MATLAB 'unique') and truncate qSUP.
%   9. Perform SVD on qSUP for nonlinear fitting (UU, SS, VV).
%
% Features & Notes
% ----------------
%   • Handles **non-homogeneous Dirichlet BCs** by consistently using the
%     unconstrained subspace (DOFl) for SVD extraction.
%   • Produces diagnostic plot of qINF vs snapshot index.
%   • Designed to feed manifold-learning or spline-regression routines
%     (e.g., Bsplines_EP_1stmodeDMG, RBF mappings).
%   • Optional tolerance control via DATAoffline.errorDISP.
%
% Reference Path
% --------------
%   /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/
%
% Author
% ------
%   Joaquín A. Hernández Ortega (JAHO)
%   CIMNE / Universitat Politècnica de Catalunya (UPC)
%
% History (versioned)
% -------------------
%   2025-11-07 — Barcelona, Spain — Comments section updated automatically by ChatGPT; no logic modified. (JAHO)
%   2025-07-27 — Borač, Serbia — Added support for non-homogeneous Dirichlet BCs. (JAHO)
%   2025-07-15 — Original implementation for homogeneous BCs. (JAHO)
%==========================================================================

%--------------------------------------------------------------------------




%  if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%      % We are interested in a global basis matrix, also including the
%      % constrained DOFs
%      SNAPdisp_ALL =cell(1,length(CASES)) ;
%  end

if nargin == 0
    load('tmp.mat')
end

% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
SNAPdisp_orthog =cell(1,length(CASES)) ;  
for iproj = 1:length(CASES)
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    if iproj == 1
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'MATPRO','MESH','OPERFE','DISP_CONDITIONS',...
            'OTHER_output') ;
        DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
        DOFl = DISP_CONDITIONS.DOFl ;
        DOFr = DISP_CONDITIONS.DOFr ;
        OTHER_output.DISP_CONDITIONS = DISP_CONDITIONS ; 
     end
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
        
        %DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        
        % Or the whole matrix ....
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        DISP_LOC{iloc} =  DISP_LOC{iloc}*SNAP_cluster.DISP.V(2:end,:)' ;
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    

    
end


% Let us begin by extracting, for instance, the x coordinates of the top
% fiber
ndim  = size(MESH.COOR,2) ; 
nDOFS  = size(MESH.COOR,2)*size(MESH.COOR,1) ; 

nodeSELECTED = DATA_interp.NodeToFollowAsLatentVariable ; 
dofLOC = DATA_interp.DirectionDisplacementAsLatentVariable ;
DOFlatent = small2large(nodeSELECTED,ndim) ;
DOFlatent = DOFlatent(dofLOC) ; 
  
 
SNAPdisp = cell2mat(SNAPdisp) ;
 
  
qMASTER=  SNAPdisp(DOFlatent,:) ; % This play the role of latent variable (master variable)
 
figure(190)
hold on
xlabel('Time step')
title(' Displacement selected DOF (latent variable), check it is injective')
ylabel('q')
plot(qMASTER)

 


% Orthogonal complement and projection onto the nonlinear basis


% 
 TOL_BLOCK = DATAoffline.errorDISP  ;
% 
% % This is just for plotting purposes. 
% % We apply the SVD to a matrix containing both unconstrained and
% % constrained DOFs
% DATAsvd=[];
% [OTHER_output.Phi_To_Plot,S,V] = SRSVD(SNAPdisp,TOL_BLOCK) ;


% SVD of the unconstrained block of SNAPdisp

      SNAPdisp =  SNAPdisp(DOFl,:) ; 
 



DATAsvd=[];
[PhiNON,S,V] = SRSVD(SNAPdisp,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes =',num2str(size(PhiNON,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])



PhiLIN = [] ; 
%PhiNON = Phi(:,2:end) ; 
% 
% % TRAIN NEURAL NETWORK, COEFFICIENTS PROJECTION
% 
%   
% 
% qSUP = cell2mat(qSUP) ;
% qSUP = qSUP(:,aaa) ;
% 

qNON = PhiNON'*SNAPdisp ;  % It may be done as well with the SVD S and V...
OTHER_output.Phi_To_Plot = zeros(nDOFS,size(PhiNON,2)) ; 
OTHER_output.Phi_To_Plot(DOFl,:) = PhiNON; 

  [UU,SS,VV ] = SVDT(qNON) ; % I know this is redundant, but we shall keep this way for compatibility with
  % previous versions


%SS = SS/SS(1) ;