function [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
    = Determine_qinf_qsup_1SVDd(CASES,NAMEsnap_base,DATAoffline)
%==========================================================================
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


% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
SNAPdisp_orthog =cell(1,length(CASES)) ;  

%  if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%      % We are interested in a global basis matrix, also including the
%      % constrained DOFs
%      SNAPdisp_ALL =cell(1,length(CASES)) ;
%  end

if nargin == 0
    load('tmp.mat')
end


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
        DISP_LOC{iloc} =  DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    

    
end

% Orthogonal complement and projection onto the nonlinear basis



TOL_BLOCK = DATAoffline.errorDISP*ones(length(SNAPdisp),1)' ;

% This is just for plotting purposes. 
% We apply the SVD to a matrix containing both unconstrained and
% constrained DOFs
DATAsvd=[];
[OTHER_output.Phi_To_Plot,S,V] = SRSVD(SNAPdisp,TOL_BLOCK,DATAsvd) ;


% SVD of the unconstrained block of SNAPdisp

for irpo = 1:length(SNAPdisp)
     SNAPdisp{irpo} =  SNAPdisp{irpo}(DOFl,:) ; 
end



DATAsvd=[];
[Phi,S,V] = SRSVD(SNAPdisp,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes =',num2str(size(Phi,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])



PhiLIN = Phi(:,1) ; 
PhiNON = Phi(:,2:end) ; 

% TRAIN NEURAL NETWORK, COEFFICIENTS PROJECTION

% Now that we know the modes BasisU  = [PhiLIN,PhiNON], we can collect the
% data for training or nonlinear fitting function (either NNs, Radial Basis functios...)
qINF = cell(size(SNAPdisp)) ;
qSUP = cell(size(SNAPdisp)) ;

for iproj = 1:length(SNAPdisp)
    qINF{iproj} = PhiLIN'*SNAPdisp{iproj} ;
    qSUP{iproj} = PhiNON'*SNAPdisp{iproj} ;
end
disp('BE CAREFUL, THIS OPERATION IS ONLY VALID IF qINF is a monotous function of input parameter')
qINF = cell2mat(qINF) ;
[qINF,aaa,bbb ]= unique(qINF) ;

figure(10)
hold on
xlabel('Number of snapshot')
ylabel('qINF')
plot(qINF,'Marker','*')
grid on

qSUP = cell2mat(qSUP) ;
qSUP = qSUP(:,aaa) ;

[UU,SS,VV ] = SVDT(qSUP) ;
%SS = SS/SS(1) ;