function [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,...
    MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,OTHER_output,...
    DATA_evaluateTAU_and_DER, nREDcoor,GmatN,DATA_interp] =Determine_qinf_qsup_DAMAG_1DK(CASES,NAMEsnap_base,DATAoffline,DATA_interp)
%==========================================================================
% File: Determine_qinf_qsup_DAMAG_1DK.m
%
% Purpose
% -------
%   End-to-end OFFLINE constructor for a 1D *damage* manifold using a
%   K_ll-based inner product. It:
%     (i)   loads and reconstructs displacement & internal-variable snapshots
%           from multiple training CASES,
%     (ii)  classifies snapshots into elastic vs nonlinear (via damage IV),
%     (iii) builds elastic and nonlinear subspaces in the G-inner product
%           (G := K_ll or M_geo,ll),
%     (iv)  selects a *nonlinear master* mode and fits B-splines for slave
%           amplitudes vs. the master coordinate,
%     (v)   assembles ROM ingredients to evaluate τ(q), τ′(q), τ″(q).
%
% Function Signature
% ------------------
%   [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl, ...
%    MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,OTHER_output, ...
%    DATA_evaluateTAU_and_DER,nREDcoor,GmatN,DATA_interp] = ...
%      Determine_qinf_qsup_DAMAG_1DK(CASES,NAMEsnap_base,DATAoffline,DATA_interp)
%
% Inputs
% ------
%   CASES         : vector of training case IDs (folder suffixes).
%   NAMEsnap_base : base path; snapshots expected at [NAMEsnap_base num2str(id)].
%   DATAoffline   : offline options (norm/metric, tolerances, mode caps, etc.).
%   DATA_interp   : B-spline learning options (order, knots, subsampling, etc.).
%
% Outputs
% -------
%   PhiMaster_lin   : [nDOF×1] elastic/linear mode (full DOFs for plotting).
%   PhiMaster_nonl  : [nDOF×1] nonlinear *master* damage mode.
%   PhiSlave_nonl   : [nDOF×(r_nl−1)] nonlinear *slave* modes.
%   MESH, DATA      : FE mesh and run metadata recovered from snapshots.
%   DOFl, DOFr      : free and constrained DOF indices.
%   OPERFE          : FE operators (quadrature, B, weights, etc.).
%   MATPRO          : material properties used in training simulations.
%   OTHER_output    : misc. info (constraint map A, K, indices, plotting bases).
%   DATA_evaluateTAU_and_DER :
%                     evaluators/metadata to compute τ(q), τ′(q), τ″(q).
%   nREDcoor        : reduced dimension of the manifold (typically 2).
%   GmatN           : inner-product matrix actually used (K_ll or M_geo,ll).
%   DATA_interp     : echoed back, possibly augmented (e.g., latent index).
%
% Workflow (high level)
% ---------------------
%   1) For each CASE:
%        • Load INFO_SNAPSHOTS and SVD blocks for displacement & internal var.
%        • Reconstruct full snapshot matrices (U*S*Vᵀ).
%   2) Concatenate all CASES → global SNAPdisp and INTV_GLO.
%   3) Classify snapshots:
%        • Elastic: mean|INTV| ≈ 0 (within tolerance).
%        • Nonlinear: remaining indices.
%   4) Call DetermineBasis_elast_damag_K to:
%        • choose G (K_ll or M_geo,ll) and build Φ_lin by WSVDT,
%        • extract nonlinear G-orthogonal complement,
%        • select nonlinear master/slave split and fit B-splines,
%        • return τ(q) evaluators and bookkeeping.
%
% Notes
% -----
%   • Internal-variable name for classification is set to 'd_DAMAGE' by default
%     (DATAoffline.NameInternalVariableSelectNonlinearSnapshots).
%   • This routine orchestrates I/O & classification; basis math and B-spline
%     fitting are delegated to DetermineBasis_elast_damag_K and its helpers.
%   • The decoder structure assumed downstream is compatible with a two-latent
%     parametrization q = (qLIN, qNON), where slave amplitudes are functions of
%     the nonlinear master coordinate.
%
% Dependencies (invoked)
% ----------------------
%   DetermineBasis_elast_damag_K, DefaultField, SVD/SRSVD utilities,
%   GiD post helpers (via OTHER_output), material/mesh readers via INFO_SNAPSHOTS.
%
% Author
% ------
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   CIMNE / Universitat Politècnica de Catalunya (UPC)
%
% History (versioned)
% -------------------
%   2025-11-07 — Barcelona, Spain — Header/comments generated automatically by ChatGPT; no logic changes. (JAHO)
%   2025-10-23 — Barcelona, Spain — Initial damage-manifold adaptation notes.
%========================================================================== 


if nargin == 0
    load('tmp1.mat')
end

% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
INTV_GLO =cell(1,length(CASES)) ;


%  if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%      % We are interested in a global basis matrix, also including the
%      % constrained DOFs
%      SNAPdisp_ALL =cell(1,length(CASES)) ;
%  end

DATAoffline.NameInternalVariableSelectNonlinearSnapshots = 'd_DAMAGE'; % 'InternalVarStrain' ;

%DATAoffline = DefaultField(DATAoffline,'NameInternalVariableSelectNonlinearSnapshots','InternalVarStrain') ;
NameINTV = DATAoffline.NameInternalVariableSelectNonlinearSnapshots ;

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
    end
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    INTV  = cell(1,length(NAME_SNAP_loc)) ; % Internal variable
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
        
        %DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        
        % Or the whole matrix ....
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        DISP_LOC{iloc} =  DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
        
        % dECOMPRESSION FROM THE SVD
        INTV{iloc} = bsxfun(@times,SNAP_cluster.(NameINTV).U',SNAP_cluster.(NameINTV).S)' ;
        INTV{iloc} =  INTV{iloc}*SNAP_cluster.(NameINTV).V' ;
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    INTV_GLO{iproj} = cell2mat(INTV);
    %     if iproj == 1
    %         % Nonlinear problems
    %         % First nonzero colum in SNAPdisp{1} is considered the linear mode
    %         isnap  =2 ;
    %         PhiLIN = SNAPdisp{iproj}(:,isnap) ;
    %         n_PhiLIN = norm(PhiLIN,'fro') ;
    %         if n_PhiLIN >0
    %             PhiLIN = PhiLIN/n_PhiLIN ;
    %         else
    %             error('Zero norm...Choose another snapshot')
    %         end
    %
    %     end
    %     SNAPdisp_orthog{iproj} = SNAPdisp{iproj} - PhiLIN*(PhiLIN'*SNAPdisp{iproj}) ;
    
end

% LET US BEGIN BY PUTTING ALL SNAPSHOTS IN A SINGLE MATRIX
SNAPdisp = cell2mat(SNAPdisp) ;  % displacement snapshots
INTV_GLO = cell2mat(INTV_GLO) ;  % Internal variable snapshots

% Now we have to discriminate between "linear/elastic" snapshots and
% nonlinear ones (remaining ones)
which_INTV = sum(abs(INTV_GLO),1)/size(INTV_GLO,1) ;
TOLLOC = 1e-12 ;
ind_linear = find(which_INTV <= TOLLOC) ;
ind_nonlinear = 1:length(which_INTV)  ;
ind_nonlinear = setdiff(ind_nonlinear,ind_linear) ;

% SVD for elastic part (unconstrained DOFs)
% ------------------------

OTHER_output.ind_nonlinear = ind_nonlinear;
OTHER_output.ind_linear = ind_linear;
OTHER_output.DISP_CONDITIONS = DISP_CONDITIONS ;

[PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
    qSLAVE_nonl,OTHER_output,DATA_evaluateTAU_and_DER, nREDcoor,DATA_interp,GmatN] =...
    DetermineBasis_elast_damag_K(SNAPdisp,DOFl,ind_linear,ind_nonlinear,...
    DATAoffline,DATA_interp,OTHER_output,OPERFE,INTV_GLO ) ;




