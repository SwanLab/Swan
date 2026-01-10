function OFFLINE_STAGE_nonECM_damagK(DATAoffline,DATA_interp,...
    NAME_BASE,DATA_GENGAUSS,InputDataForForces)
%==========================================================================
% File: OFFLINE_STAGE_nonECM_damagK.m
%
% Purpose
% -------
%   Build all OFFLINE artifacts to run a nonlinear manifold HROM for
%   small-strain isotropic DAMAGE (univariate) using a Kll-based metric and
%   ECM hyperreduction. This routine:
%     1) Extracts linear/nonlinear displacement subspaces and assembles the
%        decoder τ(q) with least-squares B-splines (per DATA_interp).
%     2) (Optionally) constructs an encoder to recover latent coordinates
%        from displacements for consistent stress evaluation.
%     3) Projects and checks PK2 stresses against a target tolerance.
%     4) Assembles internal-force snapshot matrices A_fint(qk).
%     5) Selects ECM points/weights (STANDARD or CECM-based strategies).
%     6) Stores OFFLINE data (ECMdata, BasisU, BasisStwo, DATA, …).
%
% Function Signature
% ------------------
%   OFFLINE_STAGE_nonECM_damagK(DATAoffline, DATA_interp, ...
%       NAME_BASE, DATA_GENGAUSS, InputDataForForces)
%
% Inputs
% ------
%   DATAoffline :
%     • UseEncoderToDetermineStresses (logical)
%     • errorDISP, errorSTRESS, errorPK2stress_basis
%     • errorFINT, errorFINT_linearpart
%     • nmodes_INELASTIC (cap for nonlinear modes)
%     • Hyperreduction_METHOD_projection = 'STANDARD' | 'CECM_BASED_STRATEGY'
%     • USE_ELEMENT_BASED_APPROACH (0/1), ECMforReactionForces (0/1)
%     • DATA_interp_ECM (optional override for ECM stage)
%
%   DATA_interp :
%     • METHOD_SELECT_REFERENCE_MODE (e.g., 'FIRST_SVD_MODE')
%     • METHOD_INTERP = 'BSPLINES_LEAST_SQUARES'
%     • NSAMPLES, ratio_NSAMPLES_knots, order_Bslines
%     • PortionExtrapolation_plot
%     • MAKE_SVD_AMATRIX_per_project (0/1)
%     • IndexLatentPlasticVariable (used here as latent index for damage)
%
%   NAME_BASE       : base name for training snapshot folders under ./SNAPSHOTS/
%   DATA_GENGAUSS   : options for continuous ECM (usually inactive here)
%   InputDataForForces : function handle returning training load cases
%
% Outputs (saved under ./SNAPSHOTS/NAME_BASE/OFFLINE.mat)
% -------------------------------------------------------
%   ECMdata      : ECM selection (points/elements, weights, and, for CECM,
%                  additional manifold-adaptive data)
%   BasisU       : displacement basis [Φ_master_lin, Φ_master_nonl, Φ_slave_nonl]
%   BasisStwo    : PK2-stress basis (if assembled; may be empty by design)
%   DATA         : enriched runtime data (DOFr, encoder/decoder handles, etc.)
%   DATAoffline  : snapshot/ECM/hyperreduction options used for this run
%   GmatN        : gradient-type operator used in encoder/decoder evaluation
%   BstRED_reactions : (optional) reduced operator for reaction-forces term
%
% High-Level Workflow
% -------------------
%   (A) Read training metadata and snapshots (NAME_BASE/CASE_i/…)
%   (B) Determine 1D latent manifold for DAMAGE:
%         Determine_qinf_qsup_DAMAG_1DK(...) → BasisU, τ(q), τ′(q), τ″(q)
%   (C) Reconstruct stresses; verify PK2 error ≤ DATAoffline.errorSTRESS
%   (D) Assemble A_fint(qk); (optionally) compress per project
%   (E) Compute linear-part basis (SRSVD) to stabilize ECM
%   (F) Run DiscreteECM_givenAmat (STANDARD) or CECM_based_ManifAdWeights_PLaCE
%   (G) Store OFFLINE artifacts for online HROM
%
% Notes & Recommendations
% -----------------------
%   • If PK2 error is high, tighten errorDISP, increase nmodes_INELASTIC,
%     or refine B-spline sampling/knots.
%   • MAKE_SVD_AMATRIX_per_project=1 reduces memory at the cost of an extra SVD.
%   • Element-based ECM (USE_ELEMENT_BASED_APPROACH=1) is incompatible with
%     'CECM_BASED_STRATEGY' in this implementation.
%
% Dependencies (called internally)
% --------------------------------
%   Determine_qinf_qsup_DAMAG_1DK, EncoderDamageModel_snap_v2,
%   SnapStressFromDispLOC, BasisF_from_BasisStress_PK1,
%   DiscreteECM_givenAmat, CECM_based_ManifAdWeights_PLaCE,
%   ConvertECM_points2elements, SRSVD, RSVDT, DefaultField, etc.
%
% Author
% ------
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   CIMNE / Universitat Politècnica de Catalunya (UPC)
%
% History (versioned)
% -------------------
%   2025-11-07 — Barcelona, Spain — Header/preamble generated automatically by ChatGPT; no logic changes. (JAHO)
%   2025-10-23 — Barcelona (HGs Pedralbes), Spain — Damage adaptation notes (JAHO)
%   2025-09-17 — Barcelona, Spain — CECM strategy design (internal notes)
%   2025-08-19 — Cartagena, Spain — Comments refresh in related routines
%========================================================================== 


% =========================================================================

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


%------------------------------------------------------------
delete('OFFLINE.txt')
diary 'OFFLINE.txt'
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
% SNAPSHOTS INFO
%
if  nargin == 0
    load('tmp.mat')
end
NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];



InputForces = feval(InputDataForForces) ;

CASES = 1:length(InputForces.INPUTS_PARAMETERS) ;  % Number of training projects





%   Build a 1D manifold-based reduced description that separates elastic and
%   plastic behavior from collections of displacement snapshots. The method
%   (i) identifies an elastic subspace, (ii) extracts an inelastic (plastic)
%   subspace orthogonal to the elastic one, and (iii) parameterizes the
%   plastic “slave” coordinates as a nonlinear (B‑spline) function of a
%   single plastic “master” coordinate. The result is a decoder of the form
[PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,...
    MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,OTHER_output,...
    DATA_evaluateTAU_and_DER, nREDcoor,GmatN,DATA_interp] =Determine_qinf_qsup_DAMAG_1DK(CASES,NAMEsnap_base,DATAoffline,DATA_interp) ;





BasisU = [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl] ;
%[tau,der_tau] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,zeros(2,1),DATA_evaluateTAU_and_DER) ;



NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
if  ~exist(NAME_MODES_FOLDER)
    mkdir(NAME_MODES_FOLDER)
end
NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,'DISPmodes'] ;

NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [] ;
GidPostProcessModesDOML(MESH.COOR,MESH.CN,MESH.TypeElement,OTHER_output.Phi_To_Plot,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESH.MaterialType,DATALOC) ;
OTHER_output.Phi_To_Plot = [] ;
%

%BasisU = BasisU(DOFl,:) ;




% Now we have to compute the PK2 stresses produced by these strains, and check
% that all blocks are approximated to an accuracy theshold
% DATAoffline.\errorSTRESS

STRESS_PK2_error = zeros(length(CASES),1) ;
qLATENT = cell(1,length(CASES)) ;  % Latent variables

%SNAPstressSTWOproj = cell(1,length(CASES)) ;
%SNAPstressPonePROJ = cell(1,length(CASES)) ;
A_internalFORCES_ECM = cell(1,length(CASES)) ; % For hyperreduction purposes


if ~isempty(OTHER_output.DISP_CONDITIONS.A)
    % Affine BCs, see
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
    BstRED_l = OPERFE.Bst*OTHER_output.DISP_CONDITIONS.A*BasisU ;
else
    BstRED_l = OPERFE.Bst(:,DOFl)*BasisU ;
end

% KW:Rforce --------------------------------------------------------------------------------
% Change introduced 2-Oct-2025
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
[BstRED_r,DATAoffline] = DetermineModeReactionsOUTPUT_glo(DATAoffline,MESH,OPERFE.Bst,OPERFE.wSTs,OTHER_output) ;

NORMALIZE_B_r = 1 ; 
if ~isempty(BstRED_r) && NORMALIZE_B_r == 1
    BstRED_r_normalized = BstRED_r/norm(BstRED_r,'fro');
else
    BstRED_r_normalized = [] ; 
end
% ------------------------------------------------------------------------------------------

DATAoffline = DefaultField(DATAoffline,'ECMforReactionForces',false) ; %   = true ;


DATAoffline = DefaultField(DATAoffline,'UseEncoderToDetermineStresses',1);
%DATAoffline.UseEncoderToDetermineStresses = 1;

DATA_interp = DefaultField(DATA_interp,'MAKE_SVD_AMATRIX_per_project',0) ;

switch DATAoffline.Hyperreduction_METHOD_projection
    case   'CECM_BASED_STRATEGY'
        DATA_interp.MAKE_SVD_AMATRIX_per_project = 0;
        % MAW-ECM needs to keep track of the internal work density at
        % each point of the manifold.
end
idimLAT = 1:nREDcoor ; % Latent dimensions between the set of modal coefficients
for iproj = 1:length(CASES)
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    SNAPstressSTWOproj_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SNAPstressSTWO_LOC  = cell(1,length(NAME_SNAP_loc)) ;
    
    SNAPstressPonePROJ_LOC = cell(1,length(NAME_SNAP_loc)) ;
    A_internalFORCES_ECM_LOC = cell(1,length(NAME_SNAP_loc)) ;
    
    VAR = [] ; VARint_n = [] ;
    DATA = DefaultField(DATA,'ListFieldInternalVariables',[] ) ;
    qLATENT_LOC = cell(1,length(NAME_SNAP_loc)) ;
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        % Projected displacements
        % -----------------------
        DATAoffline.UseEncoderToDetermineStresses = 1;
        if DATAoffline.UseEncoderToDetermineStresses == 0
            
            coeff = BasisU'*(Kll*SNAP_cluster.DISP.U(DOFl,:)) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            qL_extended = coeff*SNAP_cluster.DISP.V' ;
        else
            [qL_extended,qLATENT_LOC{iloc}] = EncoderDamageModel_snap_v2(BasisU(:,1:nREDcoor),GmatN,SNAP_cluster,DOFl,idimLAT,...
                DATA_evaluateTAU_and_DER) ;
            
            CHECK_eNCODER = 1; 
            if CHECK_eNCODER ==1 
                figure(3565)
                hold on
                xlabel('STEP')
                ylabel('q')
                
                CheckEncoder_here =0 ; 
                if CheckEncoder_here == 1
                 CheckEncoder(BasisU,GmatN,SNAP_cluster,qL_extended,DOFl) ;
                end
            end
            
            
            
%             
%             [qL_extended,qLATENT_LOC{iloc}] = EncoderDamageModel_snap(BasisU(:,1:nREDcoor),GmatN,SNAP_cluster,DOFl,idimLAT,...
%                 DATA_evaluateTAU_and_DER) ;
        end
        
        
        
        d = ObtainDisplacementFromSnapshotsProj(BasisU,qL_extended,SNAP_cluster,DOFr,OTHER_output,DOFl) ;
        
        %
        [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR] = SnapStressFromDispLOC...
            (VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
            SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj) ;
        
        % INTERNAL FORCES MATRIX (FOR HYPERREDUCTION PURPOSES)
        A_fint = cell(1,size(qL_extended,2)) ;
        
        for itimeLOC = 1:size(qL_extended,2)
            q =  qL_extended(1:nREDcoor,itimeLOC)  ;
            %    tauNONder_q = tauNONder(q);
            [~,tauNONder_q,~] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,q,DATA_evaluateTAU_and_DER) ;
            BstRED_l_q = BstRED_l*tauNONder_q ;
            Pk1_stress = SNAPstressPonePROJ_LOC{iloc}(:,itimeLOC);
            %  A_fint{itimeLOC} = BasisF_from_BasisStress_PK1(BstRED_l_q,Pk1_stress,DATA);
            
            
            if ~DATAoffline.ECMforReactionForces
                A_fint{itimeLOC} = BasisF_from_BasisStress_PK1(BstRED_l_q,Pk1_stress,DATA);
            else
                % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
                % 2-Oct-2025
                % KW:ReactF
                %  if isempty(DATAoffline.errorFINT_reactions)
                if ~isempty(BstRED_r_normalized)
                BstRED_r_input = BstRED_r_normalized*norm(BstRED_l_q,'fro') ;
                else
                  BstRED_r_input =   BstRED_r ; 
                end
                A_fint{itimeLOC} = BasisF_from_BasisStress_PK1([BstRED_l_q,BstRED_r_input],Pk1_stress,DATA);
                % A_react{itimeLOC} = [] ;
                % else
                %     A_fint{itimeLOC} = BasisF_from_BasisStress_PK1([BstRED_l_q],Pk1_stress,DATA);
                %     A_react{itimeLOC} = BasisF_from_BasisStress_PK1([BstRED_r],Pk1_stress,DATA);
                % end
            end
            
            
            
        end
        
        
        if DATA_interp.MAKE_SVD_AMATRIX_per_project==1
            
            A_internalFORCES_ECM_LOC{iloc} = cell2mat(A_fint) ;
        else
            A_internalFORCES_ECM_LOC{iloc} = A_fint ;
        end
        
    end
    
    
    SNAPstressSTWOproj_LOC = cell2mat(SNAPstressSTWOproj_LOC) ;   % Approximate
    SNAPstressSTWO_LOC = cell2mat(SNAPstressSTWO_LOC) ;   % exact
    
    % Check if it meets the error criterion
    STRESS_PK2_error(iproj) = norm(SNAPstressSTWO_LOC-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;
    disp(['Project = ',num2str(iproj),'; ERROR stress PK2  = ',num2str(STRESS_PK2_error(iproj))]);
    if STRESS_PK2_error(iproj) > DATAoffline.errorSTRESS
        %  dbstop('129')
        error('STRESS ERROR CRITERION NOT MET: TRY WITH A LOWER ERROR TOLERANCE FOR DISPLACEMENTS ')
    end
    
    
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
    %     RELTOL = 1e-10;  % WE only need for reconstruction
    %     [UU,SS,VV] =  SRSVD(SNAPstressSTWOproj_LOC,RELTOL) ;
    %     SNAPstressSTWOproj{iproj} = bsxfun(@times,UU',SS)' ;
    %     % PK1 stresses
    %     % ----------------
    %     SNAPstressPonePROJ_LOC = cell2mat(SNAPstressPonePROJ_LOC) ;   % exact
    %     [UU,SS,VV] =  RSVDT(SNAPstressPonePROJ_LOC) ;
    %     SNAPstressPonePROJ{iproj} = bsxfun(@times,UU',SS)' ;
    
    % internal forces for_ECM
    
    
    if DATA_interp.MAKE_SVD_AMATRIX_per_project==1
        A_internalFORCES_ECM_LOC = cell2mat(A_internalFORCES_ECM_LOC) ;
        %  [UU,SS,VV] =  RSVDT(A_internalFORCES_ECM_LOC) ;
        A_internalFORCES_ECM{iproj} = A_internalFORCES_ECM_LOC ;
        % compress the information of the internal forces
        % Note that what counts is the column space of this matrix
    else
        A_internalFORCES_ECM{iproj} = horzcat(A_internalFORCES_ECM_LOC{:}) ;
    end
    
    qLATENT{iproj} = cell2mat(qLATENT_LOC) ;
    
end

if DATA_interp.MAKE_SVD_AMATRIX_per_project == 0
    A_internalFORCES_ECM =  horzcat(A_internalFORCES_ECM{:}) ;
    
end
qLATENT = cell2mat(qLATENT) ;

% BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK2 STRESSES (FOR RECONSTRUCTION PURPOSES)
% -------------------------------------------------------------------------------------------------
% TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAPstressSTWOproj),1)' ;
%
% DATAsvd=[];
% [BasisStwo,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressSTWOproj,TOL_BLOCK,DATAsvd) ;
% disp('***********************************************************')
% disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
% disp('***********************************************************')
BasisStwo = [] ;
%
%
% % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
% % ----------------------------------------------------------------------------------------------
%
% DATAsvd=[];
% [BasisPone,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressPonePROJ,TOL_BLOCK,DATAsvd) ;
% % disp('***********************************************************')
% % disp(['Number of PK1 stress modes =',num2str(size(BasisPone,2))])
% % disp('***********************************************************')
% % The basis for the PK1 stresses is weighted with the singular values
% BasisPone = bsxfun(@times,BasisPone',S)' ;

% -------------------
% HYPERREDUCTION
% -------------------

% Determination of internal force snapshot  matrix from PK1 stresses


% if DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%
%     BasisU_all = zeros(length(DOFr)+length(DOFl),size(BasisU,2)) ;
%     BasisU_all(DOFl,:) = BasisU ;
%     BasisU_all(DOFr,:) = BasisU_r ;
%
%     BstRED_l = OPERFE.Bst*BasisU_all ;
% %
% % else
% BstRED_l = OPERFE.Bst(:,DOFl)*BasisU ;
% %
% end
DATAoffline = DefaultField(DATAoffline,'USE_ELEMENT_BASED_APPROACH',0) ;
%if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0

% ******************
% Discrete ECM
% ******************
%  if DATA_GENGAUSS.ACTIVE == 0

% Elastic basis matrix
%Ael = A_internalFORCES_ECM(:,OTHER_output.ind_linear) ;

if DATAoffline.USE_ELEMENT_BASED_APPROACH == 1
    [wST,A_internalFORCES_ECM] = ConvertECM_points2elements(DATA,A_internalFORCES_ECM,OPERFE.wSTs)  ;
else
    wST = OPERFE.wSTs;
end


DATAoffline = DefaultField(DATAoffline,'errorFINT_linearpart',1e-6) ;

% LINEAR BASIS MATRIX
TOLloc = DATAoffline.errorFINT_linearpart;
DATAsvdLOC.HIDE_OUTPUT =  1;
[Qel_nw,SSS,VVV]  = SRSVD(A_internalFORCES_ECM(:,OTHER_output.ind_linear),TOLloc,DATAsvdLOC)  ;
DATA.Matrix_Include_SVD = Qel_nw ;

switch  DATA_interp.METHOD_SELECT_REFERENCE_MODE
    case 'STANDARD_ROM'
        DATAoffline.Hyperreduction_METHOD_projection = 'STANDARD'  ;
        
end


[setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM_givenAmat(A_internalFORCES_ECM,DATA,wST,DATAoffline) ;




%   DATAoffline.Hyperreduction_METHOD_projection = 'MANIFOLD';
switch DATAoffline.Hyperreduction_METHOD_projection
    
    case 'STANDARD'
        if  DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
            ECMdata.setPoints = setPoints ;
            ECMdata.wRED = wRED ;
            proporP = length(setPoints)/DATA.MESH.ngausT*100;
            disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])
            setElements = large2smallREP(setPoints,DATA.MESH.ngaus) ;
            disp('****************************+')
            disp(['List of selected m = ',num2str(length(setElements)),' elements'])
            disp(num2str(setElements'))
            clipboard('copy',num2str(setElements'));
            ECMdata.setElements = setElements ;
        else
            disp(['Element-based ecm *************************++'])
            setElements = setPoints;
            disp(['List of selected m = ',num2str(length(setElements)),' elements'])
            disp(num2str(setElements'))
            
            setPoints = small2large(setElements,DATA.MESH.ngaus_STRESS) ;
            disp(['Total number of Gauss points = ',num2str(length(setPoints))])
            %   wRED = repmat(wRED./wST(setElements),DATA.MESH.ngaus_STRESS,1) ;
            wRED = repmat(wRED',DATA.MESH.ngaus_STRESS,1) ;
            wRED = wRED(:).*OPERFE.wSTs(setPoints,:) ;
            
            ECMdata.setPoints = setPoints ;
            ECMdata.wRED = wRED ;
            ECMdata.setElements = setElements ;
            
        end
        
    case 'CECM_BASED_STRATEGY'
        
        disp('Hyperreduction strategy based on the continuous ECM ')
        
        % Strategy devised on 17th September 2025, see
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/11_MAW_ECM_plast.mlx
        
        %                 No_enf_negat_w  =2 ;
        %                 if No_enf_negat_w == 1
        %                 ECMdata = CECM_based_ManifAdWeights_PL(A_internalFORCES_ECM,DATA,OPERFE.wSTs,DATA_interp,setPoints,wRED,qLATENT,...
        %                     OTHER_output,Qel_nw,DATAOUT.Qbasis_weighted,DATAoffline) ;
        %                 elseif No_enf_negat_w == 0
        %                     ECMdata = CECM_based_ManifAdWeights_PLa(A_internalFORCES_ECM,DATA,OPERFE.wSTs,DATA_interp,setPoints,wRED,qLATENT,...
        %                     OTHER_output,Qel_nw,DATAOUT.Qbasis_weighted,DATAoffline) ;
        %                 elseif No_enf_negat_w == 2
        DATA_interp = DefaultField(DATA_interp,'IndexLatentPlasticVariable',2) ;
        DATAoffline = DefaultField(DATAoffline,'DATA_interp_ECM',DATA_interp) ;
        DATAoffline.DATA_interp_ECM.IndexLatentPlasticVariable = DATA_interp.IndexLatentPlasticVariable ;
        
        ECMdata = CECM_based_ManifAdWeights_PLaCE(A_internalFORCES_ECM,DATA,wST,DATAoffline.DATA_interp_ECM,setPoints,wRED,qLATENT,...
            OTHER_output,Qel_nw,DATAOUT.Qbasis_weighted,DATAoffline) ;
        
        
        if  DATAoffline.USE_ELEMENT_BASED_APPROACH == 1
            error('USE_ELEMENT_BASED_APPROACH = 1 is not compatible with this option ')
        end
        
        %                 else
        %                     error('Option not implemented')
        %                 end
        
    case 'MANIFOLD'
        error('Option superseded by  CECM_BASED_STRATEGY '  )
        disp('Master/Slave nonlinear mapping between ECM points ')
        % Notice that now the "master" points play the role of the
        % actual ECM points (for purposes of allocating memory for internal variables, for instance
        % What happens at slave points in terms of stresese is not of corner for the method)
        
        DATA_interp =DefaultField(DATA_interp,'NumberOfMasterPOINTS_ECMadaptive',1);
        
        if DATA_interp.NumberOfMasterPOINTS_ECMadaptive == 1
            [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.etaNON,ECMdata.etaNONder,...
                ECMdata.etaNONder2] =...
                DiscreteECM_adaptWEIGHTS(A_internalFORCES_ECM,setPoints,wRED,DATA_interp,OPERFE.wSTs) ;
        elseif DATA_interp.NumberOfMasterPOINTS_ECMadaptive == 2
            [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.etaNON,ECMdata.etaNONder,...
                ECMdata.etaNONder2] =...
                DiscreteECM_adaptWEIGHTS_2p(A_internalFORCES_ECM,setPoints,wRED,DATA_interp,OPERFE.wSTs) ;
        else
            error('General case not implemented yet (5th August 2025)')
        end
        
        
        
        
        ECMdata.setElements = large2smallREP(ECMdata.setPoints,DATA.MESH.ngaus) ;
        disp(['Master element(s) =',num2str(ECMdata.setElements(:)')])
        ECMdata.setElements_slv = large2smallREP(ECMdata.setPoints_slv,DATA.MESH.ngaus) ;
        disp(['Slave element(s) =',num2str(ECMdata.setElements_slv(:)')])
    otherwise
        ECMdata.etaNON =[] ;
end


%
%     else
%         error('Option not available (2-July-2025)')
%         % ****************************************************************
%         % Function for computing the position of the integration points
%         % ****************************************************************
%         Nst = OTHER_output.Nst ;
%         NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
%         if  ~exist(NAME_MODES_FOLDER)
%             mkdir(NAME_MODES_FOLDER)
%         end
%         DATA_GENGAUSS.NameFileMesh_ECM = [NAME_MODES_FOLDER,NAME_BASE,'DECMpoints'] ;
%         DATALOC = [] ;
%
%         DATA_GENGAUSS.NameFileMesh_FINT = [NAME_MODES_FOLDER,NAME_BASE,'InternalForceModes'] ;
%         DATALOC = [] ;
%
%         DATA_GENGAUSS.NameFileMesh_CECM = [NAME_MODES_FOLDER,NAME_BASE,'CECMpoints'] ;
%         DATALOC = [] ;
%
%
%
%         [ECMdata] = ContinuousECM(BstRED_l,BasisPone,DATA,OPERFE.wSTs,DATAoffline,DATA_GENGAUSS,...
%             MESH,Nst) ;
%     end

% else
%     error('Option not maintained (2-Jul-2025)')
%     SNAPredFINT = BasisF_from_BasisStress_PK1_ELEMS(BstRED_l,BasisPone,DATA, OPERFE.wSTs)  ;
%     wSTs_LOC = ones(size(SNAPredFINT,1),1) ;
%     %     sqrt_wST = sqrt(OPERFE.wSTs) ;
%     %     SNAPredFINT = bsxfun(@times,SNAPredFINT,sqrt_wST) ;
%     % Determine an  orthogonal basis matrix $Q$ for the column space of $\SNAPredFINTw{}{}$
%     %DATAoffline.errorFINT = 1e-3;
%     DATAsvd.RELATIVE_SVD = 1;
%     [Q,S,V,eSVD,Rsup] = RSVDT(SNAPredFINT,DATAoffline.errorFINT,[],0,DATAsvd) ;
%
%     if DATAoffline.errorFINT == 0
%         ifig = 3000 ;
%         SVDplotERROR_local(S,ifig) ;
%     end
%
%     % % Enlarge the basis matris for SNAPredFINT
%     a  = wSTs_LOC - Q*(Q'*wSTs_LOC) ;
%     if norm(a) > 1e-10
%         a = a/norm(a) ;
%         Q = [Q,a] ;
%     end
%     % Empirical cubature method
%     % -------------------------
%     DATA_ECM = [] ;
%     DATA_ECM.TOL = DATAoffline.errorECM ;
%     [setElements,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs_LOC,DATA_ECM)  ;
%     disp(['Element-based approach *************************++'])
%     disp(['List of selected m = ',num2str(length(setElements)),' elements'])
%     disp(num2str(setElements'))
%     figure(345)
%     hold on
%     xlabel('Points')
%     ylabel('Weights')
%
%     bar(sort(wRED,'descend'))
%
%     % Determine set of points
%     setPoints = small2large(setElements,DATA.MESH.ngaus_STRESS) ;
%     disp(['Total number of Gauss points = ',num2str(length(setPoints))])
%     wRED = repmat(wRED',DATA.MESH.ngaus_STRESS,1) ;
%     wRED = wRED(:).*OPERFE.wSTs(setPoints,:) ;
%
%     ECMdata.setPoints = setPoints ;
%     ECMdata.wRED = wRED ;
%     ECMdata.setElements = setElements ;

%end

disp(['*********************************************************************'])
disp(['Number of displacement modes = ',num2str(size(BasisU,2))])
disp(['Number of PK2-stress modes = ',num2str(size(BasisStwo,2))])
%disp(['Number of PK1-stress modes = ',num2str(size(BasisPone,2))])



% STORING INFORMATION ---HYPERREDUCED-ORDER OPERATORS

%BASES.BasisU = BasisU ;
DATA.DOFr = DOFr;
DATA.DATA_evaluateTAU_and_DER = DATA_evaluateTAU_and_DER;
if ~isempty(BstRED_r)
    setEntries = small2large(ECMdata.setPoints,DATA.MESH.nstrain_F) ;
    BstRED_reactions = BstRED_r(setEntries,:) ;
else
    BstRED_reactions  = [] ;
end
DATA.nREDcoor = nREDcoor;


% switch  DATAoffline.Hyperreduction_METHOD_projection
%     case 'STANDARD'
%         DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0;
%
%     case  'STANDARD_NONLINEAR_PART_mixed'
%         DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 1;
%
%
% end

%KlinLLred = BasisU'*(Kll*BasisU) ;
%Klin = OTHER_output.K;  % This is redundant...In future versions, do not store Kll = Klin(DOFl,
save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','GmatN','BstRED_reactions' )




%save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','Kll')



diary off

