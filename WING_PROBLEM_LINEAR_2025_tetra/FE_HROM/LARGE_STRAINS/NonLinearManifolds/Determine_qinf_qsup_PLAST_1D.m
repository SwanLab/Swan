function [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,...
    MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,OTHER_output,...
    DATA_evaluateTAU_and_DER, nREDcoor] =Determine_qinf_qsup_PLAST_1D(CASES,NAMEsnap_base,DATAoffline,DATA_interp)
%--------------------------------------------------------------------------
% function [PhiMaster_lin, PhiMaster_nonl, PhiSlave_nonl, ...
%           MESH, DATA, DOFl, DOFr, OPERFE, MATPRO, OTHER_output, ...
%           DATA_evaluateTAU_and_DER, nREDcoor] = ...
%           Determine_qinf_qsup_PLAST_1D(CASES, NAMEsnap_base, DATAoffline, DATA_interp)
%
% PURPOSE
%   Constructs a 1D manifold-based reduced representation for elastoplastic
%   problems. The method:
%     1) Separates elastic and plastic contributions from displacement snapshots
%        using the internal variable (e.g., plastic strain) as an indicator.
%     2) Builds an elastic reduced basis (PhiMaster_lin).
%     3) Identifies a plastic subspace orthogonal to the elastic one, with
%        the leading mode as the "master" plastic coordinate (PhiMaster_nonl)
%        and the remaining modes as "slave" plastic coordinates (PhiSlave_nonl).
%     4) Learns a nonlinear mapping qSLAVE_nonl = f(qMASTER_nonl) using
%        least-squares B-splines in a compact SVD basis.
%
%   The final reduced representation is:
%
%       d_L ≈ PhiMaster_lin * qMASTER_lin
%            + PhiMaster_nonl * qMASTER_nonl
%            + PhiSlave_nonl  * f(qMASTER_nonl),
%
%   where f(·) is provided by spline-based function handles and derivatives.
%
% INPUTS
%   CASES         : Vector of case IDs (each corresponds to a snapshot folder).
%   NAMEsnap_base : Base string for snapshot directories (e.g. 'DATA_snap_case_').
%   DATAoffline   : Struct with offline tolerances and options (e.g., errorDISP).
%   DATA_interp   : Struct with spline fitting options. Typical fields:
%                     • NSAMPLES                   - number of B-spline breaks
%                     • INCLUDE_SECOND_DERIVATIVES - whether to fit f'' as well
%                     • EXTRAPOLATION              - strategy outside training range
%                     • Orthogonality_elast_plast  - 'Euclidean' (default) or 'K'
%
% OUTPUTS
%   PhiMaster_lin       : [nDOF_L × r_el] basis for elastic subspace (left DOFs only).
%   PhiMaster_nonl      : [nDOF_L × 1] master plastic mode, orthogonal to PhiMaster_lin.
%   PhiSlave_nonl       : [nDOF_L × (r_pl-1)] slave plastic modes.
%   MESH                : Mesh structure loaded from snapshots.
%   DATA                : Simulation data used in snapshots.
%   DOFl, DOFr          : Index sets of left/right constrained DOFs.
%   OPERFE              : Finite element operators (assembly, strain-displacement, etc.).
%   MATPRO              : Material property structure (includes plasticity data).
%   OTHER_output        : Additional auxiliary data stored with snapshots.
%   DATA_evaluateTAU_and_DER : Struct containing function handles to evaluate
%                              f(qMASTER_nonl), f'(qMASTER_nonl), and optionally f''.
%   nREDcoor            : Effective reduced dimensionality for the plastic mapping.
%
% ALGORITHM
%   1) Load displacement and internal variable snapshots for all CASES.
%   2) Concatenate them into global matrices.
%   3) Identify elastic vs plastic snapshots using internal variable threshold (TOLLOC).
%   4) Perform SVD on elastic snapshots (DOFl only) → PhiMaster_lin.
%   5) Project plastic snapshots orthogonal to PhiMaster_lin and compute SVD:
%        - First plastic mode → PhiMaster_nonl
%        - Remaining plastic modes → PhiSlave_nonl
%   6) Project snapshots onto plastic basis → (qMASTER_nonl, qSLAVE_nonl).
%   7) Sort by qMASTER_nonl (monotone abscissa).
%   8) Compute compact SVD of qSLAVE_nonl and fit splines to right singular
%      vectors as functions of qMASTER_nonl.
%   9) Return spline handles (f, f', f'') and reduced dimensionality nREDcoor.
%
% NOTES
%   • Elastic/plastic discrimination is done via mean absolute value of
%     internal variable per snapshot.
%   • Sorting by qMASTER_nonl before spline fitting is essential.
%   • Compact SVD reduces number of spline fits and improves conditioning.
%   • DATA_interp.Orthogonality_elast_plast = 'K' can enforce stiffness-based
%     inner product instead of Euclidean.
%
% EXAMPLE
%   DATAoffline.errorDISP = 1e-10;
%   DATA_interp.NSAMPLES = 20;
%   [PhiMaster_lin, PhiMaster_nonl, PhiSlave_nonl, MESH, DATA, DOFl, DOFr, ...
%    OPERFE, MATPRO, OTHER_output, DATA_evaluateTAU_and_DER, nREDcoor] = ...
%        Determine_qinf_qsup_PLAST_1D(CASES, NAMEsnap_base, DATAoffline, DATA_interp);
%
% VERSION HISTORY
%   • 2025-07-27 (Borac, Serbia) - Plastic/elastic split version with 1D manifold mapping (J.A. Hernández).
%   • Based on Determine_qinf_qsup_1SVD, adapted August 2025 (Molinos Marfagones, Cartagena).
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
        INTV{iloc} = bsxfun(@times,SNAP_cluster.InternalVarStrain.U',SNAP_cluster.InternalVarStrain.S)' ;
        INTV{iloc} =  INTV{iloc}*SNAP_cluster.InternalVarStrain.V' ;
        
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


DATA_interp = DefaultField(DATA_interp,'Orthogonality_elast_plast','Euclidean') ; % = 'K' ;

switch DATA_interp.Orthogonality_elast_plast
    case 'Euclidean'
        [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
            qSLAVE_nonl,OTHER_output] = DetermineBasis_elast_plast_EUCLID(SNAPdisp,DOFl,ind_elastic,ind_plastic,...
            DATAoffline,DATA_interp,OTHER_output) ;
        
        
        % Since we intend to construct qSLAVE_nonl = f(qMASTER_nonl) using
        % B-splines, it is necessary to sort the entries of qMASTER_nonl is
        % ascending order
        [qMASTER_nonl,indORDER] = sort(qMASTER_nonl) ;
        
        figure(126)
        hold on
        xlabel('Snapshot ordered')
        ylabel('Master function (qMASTER nonl)')
        plot(qMASTER_nonl)
        
        qSLAVE_nonl = qSLAVE_nonl(:,indORDER)  ;
        
        % Now we have to establish the mapping qSLAVE_nonl = f(qMASTER_nonl)
        % We do this in an indirect way, by first making
        %% qSLAVE_nonl = UU*diag(SS)*VV^T
        % and then constructing the mapping VV(:,i) = f(qMASTER_nonl)
        %
        %
        
        %
        
        
        [UU,SS,VV] = SVDT(qSLAVE_nonl) ;
        
        % figure(345)
        % hold on
        % title('Right singular vectors qSLAVE nonl as a function of qMASTER nonl')
        % xlabel('qMASTER nonl')
        % ylabel('qSLAVE nonl')
        %
        %
        % for iii = 1:size(VV,2)
        %     plot(qMASTER_nonl,VV(:,iii),'DisplayName',['\lambda',num2str(iii),'=',num2str(SS(iii)/SS(1))]) ;
        % end
        
        
        [DATA_evaluateTAU_and_DER, nREDcoor] = BsplinesLeastSquares_PLAST(DATA_interp, qMASTER_nonl, VV', UU, SS) ;
        
        
        
        
    case 'K'
        
        [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
            qSLAVE_nonl,OTHER_output,DATA_evaluateTAU_and_DER, nREDcoor] = DetermineBasis_elast_plast_K(SNAPdisp,DOFl,ind_elastic,ind_plastic,...
            DATAoffline,DATA_interp,OTHER_output,OPERFE,INTV_GLO) ;
        
end



