function [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,...
    MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,OTHER_output,...
    DATA_evaluateTAU_and_DER, nREDcoor,Kll,DATA_interp] =Determine_qinf_qsup_P1D_DD(CASES,NAMEsnap_base,DATAoffline,DATA_interp)
% Determine_qinf_qsup_P1D_DD is a modification of
% Determine_qinf_qsup_PLAST_1DK, described in detail below. 
% The modification is intented to discover automatically the latent space
% of two variables (elastic/plastic)
%  See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/113_DDrivenDiscov_INTV/01_plastUNI.mlx
%  JAHO, 2-Nov-2025, Sunday, Balmes 185, Barcelona

% =========================================================================
% DETERMINE_QINF_QSUP_PLAST_1DK — Plastic Manifold Reduction with Kll Metric
% =========================================================================
% PURPOSE
%   Constructs a 1D manifold-based reduced representation for small-strain
%   elastoplasticity, using displacement and internal-variable snapshots
%   from multiple training cases. The method separates elastic and plastic
%   subspaces, then parameterizes plastic evolution with a *single master
%   coordinate* (first plastic SVD mode). Slave plastic coordinates are
%   expressed as nonlinear spline functions of this master variable.
%
%   The Kll stiffness matrix (constrained stiffness) is used to define all
%   inner products. This ensures energy-based orthogonality and stabilizes
%   the elastic–plastic decomposition.
%
% REDUCED REPRESENTATION
%   d_L(q) ≈ PhiMaster_lin * q_ELAST
%          + PhiMaster_nonl * q_PLAST
%          + PhiSlave_nonl  * f(q_PLAST),
%
%   where:
%     • q_ELAST  : elastic generalized coordinate,
%     • q_PLAST  : master plastic generalized coordinate (first plastic mode),
%     • f(·)     : nonlinear spline mapping providing plastic “slave” coords.
%
% INPUTS
%   CASES         : Vector of integer case IDs (snapshot sets to process).
%   NAMEsnap_base : Base directory string for snapshot folders.
%   DATAoffline   : Struct of offline tolerances/options:
%                     • errorDISP   – tolerance for displacement SVD.
%                     • nmodes_PLASTIC – max number of plastic modes kept.
%                     • errorSTRESS – reconstruction tolerance check.
%   DATA_interp   : Struct with interpolation options for f(q):
%                     • METHOD_INTERP = 'BSPLINES_LEAST_SQUARES' (default)
%                     • NSAMPLES, order_Bsplines, INCLUDE_SECOND_DERIVATIVES
%                     • Orthogonality_elast_plast = 'K' (enforce Kll metric)
%
% OUTPUTS
%   PhiMaster_lin : Basis for elastic subspace (DOFl only).
%   PhiMaster_nonl: First plastic mode = master plastic coordinate.
%   PhiSlave_nonl : Remaining plastic modes (slaves).
%   MESH          : Mesh structure recovered from snapshots.
%   DATA          : FE data struct (stores material, BCs, etc.).
%   DOFl, DOFr    : Index sets for unconstrained (left) and constrained (right) DOFs.
%   OPERFE        : FE operators (B-matrix, quadrature, etc.).
%   MATPRO        : Material property struct (plasticity law).
%   OTHER_output  : Extra fields from snapshot storage (e.g. plotting modes).
%   DATA_evaluateTAU_and_DER : Handles to decoder τ(q), τ′(q), τ″(q).
%   nREDcoor      : Reduced dimensionality of plastic manifold.
%   Kll           : Stiffness matrix restricted to DOFl, used for inner products.
%
% ALGORITHM
%   1. Load displacement and internal variable (plastic strain) snapshots
%      for each CASE. Store as matrices.
%   2. Concatenate across cases into global matrices.
%   3. Classify snapshots into elastic vs plastic:
%         – Elastic: mean internal variable ≈ 0
%         – Plastic: all others
%   4. Compute SVD of elastic snapshots → PhiMaster_lin.
%   5. Project plastic snapshots orthogonal to PhiMaster_lin (using Kll).
%   6. Compute SVD of projected plastic snapshots:
%         – First mode → PhiMaster_nonl (chosen master plastic variable),
%         – Remaining → PhiSlave_nonl.
%   7. Project snapshots onto plastic basis → q_MASTER, q_SLAVE.
%   8. Sort by q_MASTER (monotone abscissa).
%   9. Apply compact SVD to q_SLAVE, fit B-splines vs q_MASTER.
%  10. Construct handles f(q), f′(q), f″(q) for decoder + derivatives.
%
% NOTES
%   • The use of the *first plastic SVD mode* as master is mandatory: other
%     candidate latent variables (e.g. force, arc-length) cannot define an
%     invertible encoder, preventing stress recovery and hyperreduction.
%   • Sorting by q_MASTER is essential to avoid non-monotone fits.
%   • Compact SVD minimizes the number of spline fits while retaining accuracy.
%   • Using the Kll metric for inner products provides physical orthogonality
%     (strain-energy sense), improving stress reconstruction robustness.
%
% VERSION HISTORY
%   • Based on Determine_qinf_qsup_PLAST_1D (27-JUL-2025, Borac, Serbia).
%   • Modified to use Kll inner products and enforced master mode choice.
%     J.A. Hernández Ortega — 17-AUG-2025, Santa Lucía Hospital, Cartagena.
% =========================================================================

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
INTV_GLO =cell(1,length(CASES)) ;
PLstrain_GLO =cell(1,length(CASES)) ;



 

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
    PLstrain_LOC  = cell(1,length(NAME_SNAP_loc)) ; % Internal variable
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
        
         PLstrain_LOC{iloc} = bsxfun(@times,SNAP_cluster.PlasticStrains.U',SNAP_cluster.PlasticStrains.S)' ;
        PLstrain_LOC{iloc} =  PLstrain_LOC{iloc}*SNAP_cluster.PlasticStrains.V' ;
        
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    INTV_GLO{iproj} = cell2mat(INTV);
        PLstrain_GLO{iproj} = cell2mat(PLstrain_LOC);

    
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
PLstrain_GLO = cell2mat(PLstrain_GLO) ;  % Internal variable snapshots


% Normalization 
SNAPdisp = SNAPdisp/norm(SNAPdisp,'fro') ; 
PLstrain_GLO = PLstrain_GLO/norm(PLstrain_GLO,'fro') ; 
% Augmented matrix 
MatrixSnapAUG =     [SNAPdisp; PLstrain_GLO] ; 
[PhiAugm,SSS,Vaugm] = SRSVD(MatrixSnapAUG,DATAoffline.errorDISP) ;





figure(666)
hold on
xlabel('Step')
ylabel('Right singular vectors')
title('Right Singular Vectors (modal amplitude augmented modes DISP + PLASTIC STRAINS)')
for iii= 1:size(Vaugm,2)
    plot(Vaugm(:,iii),'DisplayName',['q',num2str(iii)])
end
legend show


n = size(Vaugm,2) ; 

[y,c] = ChooseMonotonicLinearCombination_n(Vaugm,n)

Vmon = Vaugm*c;  


figure(667)
hold on
xlabel('Step')
ylabel('Right singular vectors')
title('Most monotonic linear combination of RSV')
for iii= 1:size(Vmon,2)
    plot(Vmon(:,iii),'DisplayName',['q',num2str(iii)])
end
legend show



error('Developed until here')


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




