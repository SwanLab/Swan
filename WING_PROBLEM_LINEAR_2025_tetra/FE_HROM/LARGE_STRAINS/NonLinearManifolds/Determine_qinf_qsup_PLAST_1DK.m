function [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,...
    MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,OTHER_output,...
    DATA_evaluateTAU_and_DER, nREDcoor,Kll,DATA_interp] =Determine_qinf_qsup_PLAST_1DK(CASES,NAMEsnap_base,DATAoffline,DATA_interp)
% =========================================================================
% DETERMINE_QINF_QSUP_PLAST_1DK — Plastic Manifold Reduction with Kll Metric
% =========================================================================
% PURPOSE
%   Build a 1D manifold-based reduced description for small-strain
%   elastoplasticity from multi-case snapshots (displacements + internal
%   variables). The procedure:
%     1) Splits snapshots into ELASTIC vs PLASTIC using an internal-variable
%        indicator; computes an elastic basis with Kll-orthogonality.
%     2) Projects plastic snapshots orthogonal to the elastic subspace
%        (Kll metric), then performs SVD:
%           – 1st plastic SVD mode = MASTER plastic coordinate q_PLAST,
%           – Remaining plastic modes = SLAVE coordinates.
%     3) Fits B-spline maps for SLAVE coordinates as functions of q_PLAST:
%           τ(q_PLAST), τ′(q_PLAST), τ″(q_PLAST).
%
% REDUCED REPRESENTATION
%   d_L(q) ≈ Φ_ELAST*q_ELAST
%           + Φ_PLAST_master*q_PLAST
%           + Φ_PLAST_slave * f(q_PLAST),
%   where f(·) is a least-squares B-spline fit.
%
% WHY Kll?
%   Inner products use the constrained stiffness Kll (strain-energy metric),
%   which enforces physically meaningful orthogonality and stabilizes the
%   elastic–plastic split and stress recovery.
%
% INPUTS
%   CASES         : Row/col vector of integer case IDs to process.
%   NAMEsnap_base : Base path to snapshot folders (e.g. ./SNAPSHOTS/NAME_BASE).
%   DATAoffline   : Offline options/tolerances (selected fields):
%                     • errorDISP        – SVD tolerance for displacement basis.
%                     • nmodes_PLASTIC   – cap on number of plastic modes.
%                     • errorSTRESS      – PK2 reconstruction check (relative).
%                     • NameInternalVariableSelectNonlinearSnapshots
%                                            ('InternalVarStrain' by default).
%   DATA_interp   : Interpolation options for plastic mapping τ(·):
%                     • METHOD_INTERP = 'BSPLINES_LEAST_SQUARES'
%                     • NSAMPLES, ratio_NSAMPLES_knots, order_Bsplines
%                     • INCLUDE_SECOND_DERIVATIVES (logical)
%                     • Orthogonality_elast_plast = 'K' (use Kll)
%
% OUTPUTS
%   PhiMaster_lin    : Elastic basis (on DOFl).
%   PhiMaster_nonl   : Master plastic mode (1st plastic SVD mode).
%   PhiSlave_nonl    : Plastic slave modes.
%   MESH             : Mesh structure from snapshots.
%   DATA             : FE data (materials, BCs, quadrature, etc.).
%   DOFl, DOFr       : Indices of free (left) and constrained (right) DOFs.
%   OPERFE           : FE operators (B-matrix, quadrature weights, etc.).
%   MATPRO           : Material properties (J2 plasticity).
%   OTHER_output     : Extra info (indices elastic/plastic, plotting helpers…).
%   DATA_evaluateTAU_and_DER : Handles to τ(q), τ′(q), τ″(q).
%   nREDcoor         : Number of reduced coordinates used downstream.
%   Kll              : Constrained stiffness metric on DOFl (inner products).
%   DATA_interp      : Possibly enriched interpolation settings (returned).
%
% ALGORITHM (high level)
%   1) Load & decompress per-case SVD snapshot blocks (U*S*V^T) for
%      displacements and internal variables.
%   2) Concatenate across cases → global snapshot matrices.
%   3) Classify elastic/plastic via mean |internal variable|.
%   4) Kll-SVD on elastic snapshots → Φ_ELAST.
%   5) Kll-orthogonalize plastic snapshots to Φ_ELAST.
%   6) SVD on plastic residual → master/slave plastic modes.
%   7) Project snapshots → q_MASTER, q_SLAVE; sort by q_MASTER (monotone abscissa).
%   8) (Optional) Compact SVD on q_SLAVE; fit B-splines vs q_MASTER.
%   9) Build callable τ, τ′, τ″ for the decoder.
%
% PRACTICAL NOTES
%   • The master plastic variable MUST be the first plastic SVD mode to
%     ensure an invertible encoder for consistent stress evaluation and
%     hyperreduction downstream.
%   • Sorting by q_MASTER is crucial to obtain stable, non-oscillatory
%     spline fits.
%   • Tighten errorDISP or increase nmodes_PLASTIC if PK2 reconstruction
%     fails later in the OFFLINE checks.
%
% DEPENDENCIES (invoked here or downstream)
%   DetermineBasis_elast_plast_K, DefaultField
%   (downstream consumers: SnapStressFromDispLOC, DiscreteECM_givenAmat, etc.)
%
% VERSION / AUTHORSHIP
%   • Based on Determine_qinf_qsup_PLAST_1D — 27-JUL-2025, Borac (Serbia).
%   • Modified to use Kll inner products & enforced 1st plastic SVD as master —
%     17-AUG-2025, Santa Lucía Hospital (Cartagena).
%   • Comments refreshed and clarified — 07-NOV-2025, Barcelona.
%
%   Author: Joaquín Alberto Hernández Ortega (JAHO) 
%
%   This header was generated automatically by ChatGPT on 07-NOV-2025.
% =========================================================================

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

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

     %   DATAoffline.NameInternalVariableSelectNonlinearSnapshots = 'd_DAMAGE';  'InternalVarStrain' ; 

DATAoffline = DefaultField(DATAoffline,'NameInternalVariableSelectNonlinearSnapshots','InternalVarStrain') ; 
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
ind_elastic = find(which_INTV <= TOLLOC) ;
ind_plastic = 1:length(which_INTV)  ;
ind_plastic = setdiff(ind_plastic,ind_elastic) ;

% SVD for elastic part (unconstrained DOFs)
% ------------------------

OTHER_output.ind_plastic = ind_plastic; 
OTHER_output.ind_elastic = ind_elastic; 
OTHER_output.DISP_CONDITIONS = DISP_CONDITIONS ; 

[PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
    qSLAVE_nonl,OTHER_output,DATA_evaluateTAU_and_DER, nREDcoor,Kll,DATA_interp] = DetermineBasis_elast_plast_K(SNAPdisp,DOFl,ind_elastic,ind_plastic,...
    DATAoffline,DATA_interp,OTHER_output,OPERFE,INTV_GLO ) ;




