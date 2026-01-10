function [PhiLIN,PhiNON,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
    = Determine_qinf_qsup_STDROM(CASES,NAMEsnap_base,DATAoffline)
%==========================================================================
% File: Determine_qinf_qsup_STDROM.m
%
% Purpose
% -------
%   Build a standard (linear) ROM basis from multiple displacement
%   snapshot sets and return the mode matrix Φ for the unconstrained DOFs.
%   Unlike the manifold versions, **all modal coefficients are treated as
%   independent generalized coordinates (q)** — there is no nonlinear
%   decoder/encoder split. Constrained DOFs are handled by working on DOFl.
%
% Context
% -------
%   This routine is a streamlined adaptation of `Determine_qinf_qsup_1SVDd`
%   for **standard ROMs** where the latent space equals the linear modal
%   subspace:
%
%       u(x,t) ≈ Φ · q(t) ,   with  Φ = PhiLIN  and  PhiNON = [] .
%
%   It aggregates snapshots across training cases, reconstructs U*S*Vᵀ
%   blocks, restricts to DOFl (non-homogeneous Dirichlet ready), and
%   computes Φ via (randomized) SVD with user-given tolerance.
%
% Inputs
% ------
%   CASES         : Vector of integer IDs of training folders to process.
%   NAMEsnap_base : Base path; case i is under [NAMEsnap_base num2str(CASES(i))].
%   DATAoffline   : Struct of offline options, notably:
%                     • errorDISP — SVD truncation tolerance for displacements.
%
% Outputs
% -------
%   PhiLIN          : [nDOFl × r] ROM basis spanning the linear subspace.
%   PhiNON          : [] (empty; not used in the standard ROM path).
%   MESH, DATA      : Mesh and metadata (from first case).
%   DOFl, DOFr      : Indices of free / constrained DOFs.
%   OPERFE          : FE operators (B, quadrature, etc.).
%   MATPRO          : Material model data used in snapshots.
%   DISP_CONDITIONS : Dirichlet info (stored for downstream assembly).
%   OTHER_output    : Aux data (e.g., Phi_To_Plot for full DOFs visualization).
%
% Method (high level)
% -------------------
%   1) For each CASE:
%        • Load INFO_SNAPSHOTS, reconstruct displacement blocks as U*S*Vᵀ.
%   2) Concatenate all reconstructed blocks into a global list.
%   3) Compute a visualization basis on full DOFs (OTHER_output.Phi_To_Plot).
%   4) Restrict each block to DOFl and compute Φ = SRSVD(SNAPdisp_DOFl, errorDISP).
%   5) Set PhiLIN = Φ, PhiNON = [].
%
% Notes
% -----
%   • Non-homogeneous Dirichlet BCs are supported by restricting to DOFl.
%   • The singular values are printed for traceability and quality control.
%   • This routine is intended for **projection-only** ROMs; hyperreduction
%     (ECM/CECM) is configured elsewhere.
%
% Provenance / Paths
% ------------------
%   Snapshot structure and storage follow the same conventions as
%   `Determine_qinf_qsup_1SVDd` (see FE_HROM testing problems under:
%   /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/)
%
% Author
% ------
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% History (versioned)
% -------------------
%   2025-11-07 — Barcelona, Spain — Comments updated automatically by ChatGPT; no code changes. (JAHO)
%   2025-10-09 — Barcelona (HGs, Pedralbes) — Initial STDROM variant created. (JAHO)
%   2025-07-27 — Borač, Serbia — Non-homogeneous Dirichlet support added in base routine. (JAHO)
%==========================================================================


if nargin == 0
    load('tmp1.mat')
end
% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
 
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



PhiLIN = Phi;  ; 
PhiNON = [] ; % Phi(:,2:end) ; 

% TRAIN NEURAL NETWORK, COEFFICIENTS PROJECTION

% Now that we know the modes BasisU  = [PhiLIN,PhiNON], we can collect the
% data for training or nonlinear fitting function (either NNs, Radial Basis functios...)
% qINF = cell(size(SNAPdisp)) ;
% qSUP = cell(size(SNAPdisp)) ;
% 
% for iproj = 1:length(SNAPdisp)
%     qINF{iproj} = PhiLIN'*SNAPdisp{iproj} ;
%     qSUP{iproj} = PhiNON'*SNAPdisp{iproj} ;
% end
% 
% qINF = cell2mat(qINF) ;
% [qINF,aaa,bbb ]= unique(qINF) ;

% figure(10)
% hold on
% xlabel('Number of snapshot')
% ylabel('qINF')
% plot(qINF,'Marker','*')
% % grid on
% 
% qSUP = cell2mat(qSUP) ;
% qSUP = qSUP(:,aaa) ;
% 
% [UU,SS,VV ] = SVDT(qSUP) ;
%SS = SS/SS(1) ;