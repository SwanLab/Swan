function OFFLINE_STAGE_nonECM_plastDD(DATAoffline,DATA_interp,...
    NAME_BASE,DATA_GENGAUSS,InputDataForForces)
%OFFLINE_STAGE_nonECM_plastDD is a modification of
%OFFLINE_STAGE_nonECM_plastK, described below.
% Its goal is to automatically discovert the latent space (that should be two variables, one elastic
% the other plastic)
% JAHO, 2-Nov-2025, Sunday. Balmes 185, Barcelona 

% =========================================================================
% OFFLINE_STAGE_nonECM_plastK: OFFLINE STAGE — NONLINEAR MANIFOLD HROM (PLASTICITY) WITH Kll METRIC
% =========================================================================
% PURPOSE
%   Build all OFFLINE artifacts needed to run a nonlinear manifold HROM for
%   small‑strain J2 plasticity on the plate‑with‑hole benchmark. The routine:
%     1) Extracts elastic and plastic displacement bases using the Kll
%        (constrained stiffness) inner product for physically meaningful
%        orthogonality.
%     2) Chooses the FIRST PLASTIC SVD MODE as the master latent coordinate
%        (q_PLAST) and fits B‑spline mappings for the plastic slave
%        coordinates → decoder τ(q).
%     3) Constructs an ENCODER (displacements → latent coordinates), required
%        to recover stresses and assemble hyperreduction operators.
%     4) Projects stresses and builds internal‑force snapshot matrices.
%     5) Selects Empirical Cubature (ECM) points/weights for hyperreduction.
%
% CONTEXT / RATIONALE
%   • Using Kll to define inner products stabilizes the elastic/plastic split
%     (strain‑energy metric) and improves stress accuracy.
%   • The master latent variable MUST be the first plastic SVD mode so that an
%     encoder exists; options such as “linear in force” cannot be inverted
%     consistently and break stress recovery.
%   • Decoder form used downstream:
%         d_L(q) = Φ_ELAST*q_ELAST +  Φ_MaSTER_plas*q_PLAST + Φ_slave_plast * f(q_PLAST)
%     where f(·) is fitted with least‑squares B‑splines (order set in
%     DATA_interp).
%
% INPUTS
%   DATAoffline : struct controlling truncations and checks, e.g.
%       • UseEncoderToDetermineStresses (0/1)   : use encoder for stress eval.
%       • errorDISP     (scalar)                : SVD tol for displacement basis.
%       • nmodes_PLASTIC (int)                  : cap on # plastic modes.
%       • errorFINT     (scalar)                : ECM tolerance for internal forces.
%       • errorSTRESS   (scalar)                : PK2 stress reconstr. check.
%       • errorPK2stress_basis (scalar)         : tol for PK2 stress basis build.
%       • Hyperreduction_METHOD_projection      : 'STANDARD' or 'MANIFOLD'.
%
%   DATA_interp : struct for manifold fitting and ECM assembly, e.g.
%       • METHOD_SELECT_REFERENCE_MODE = 'FIRST_SVD_MODE'
%       • METHOD_INTERP = 'BSPLINES_LEAST_SQUARES'
%       • NSAMPLES, ratio_NSAMPLES_knots, order_Bslines, PortionExtrapolation_plot
%       • MAKE_SVD_AMATRIX_per_project (0/1) : compress A_fint per project.
%       • LocalToleranceECMnonlinear_masterPOINTS : tol for manifold‑ECM option.
%
%   NAME_BASE   : base name for locating training snapshots under ./SNAPSHOTS.
%   DATA_GENGAUSS : (optional) options for continuous ECM (usually inactive here).
%   InputDataForForces : function handle that returns the training load cases.
%
% OUTPUTS (SAVED to ./SNAPSHOTS/NAME_BASE/OFFLINE.mat)
%   ECMdata      : selected ECM Gauss points, weights, and (if 'MANIFOLD')
%                  the master/slave mapping η,η′,η″ between ECM points.
%   BasisU       : displacement basis [Φ_ELAST, Φ_PLAST_master, Φ_PLAST_slave].
%   BasisStwo    : PK2 stress basis for reconstruction/verification.
%   DATA         : enriched runtime data (DOFr, encoder/decoder handles, etc.).
%   DATA.DATA_evaluateTAU_and_DER : callable to evaluate τ(q), τ′(q), τ″(q).
%   DATA.nREDcoor: number of generalized coordinates (q_ELAST, q_PLAST, ...).
%   Kll          : constrained stiffness matrix used as inner‑product metric.
%
% HIGH‑LEVEL WORKFLOW
%   (A) Load training cases and snapshot metadata.
%   (B) Determine_qinf_qsup_PLAST_1DK(...)
%         → Kll‑orthogonal split into elastic/plastic subspaces,
%           first plastic SVD mode = master, B‑spline decoder + encoder.
%   (C) Plot/save modes for inspection (GiD).
%   (D) Using the encoder/decoder:
%         • Reconstruct displacements (d_L,d_R) per snapshot,
%         • Compute PK2/PK1 stresses via constitutive routine,
%         • Assemble reduced internal‑force snapshots A_fint(q_k).
%   (E) Build stress and A_fint bases (RSVD/blocked SVD).
%   (F) Run Discrete ECM on A_fint to select Gauss points + weights.
%   (G) (Optional) If 'MANIFOLD' hyperreduction: compute η, η′, η″ between
%       master/slave ECM points to further compress integration.
%   (H) Sanity checks (PK2 reconstruction error ≤ errorSTRESS), store OFFLINE.
%
% PRACTICAL NOTES
%   • If PK2 stress error exceeds DATAoffline.errorSTRESS, tighten errorDISP,
%     increase nmodes_PLASTIC, or refine NSAMPLES/knots in DATA_interp.
%   • MAKE_SVD_AMATRIX_per_project = 1 greatly reduces memory for ECM
%     assembly at the cost of an extra SVD step.
%   • The encoder is mandatory for consistent stress evaluation and
%     hyperreduction. Do not disable unless debugging.
%
% DEPENDENCIES (called internally)
%   Determine_qinf_qsup_PLAST_1DK, GidPostProcessModesDOML,
%   SnapStressFromDispLOC, RSVDqp/RSVDT/SRSVD, DiscreteECM_givenAmat,
%   DiscreteECM_adaptWEIGHTS(_2p), BasisF_from_BasisStress_PK1, DefaultField
%
% VERSION / AUTHOR
%   JAHO — Comments updated 19‑AUG‑2025 (Cartagena). Replaces older headers.
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
    DATA_evaluateTAU_and_DER, nREDcoor,Kll,DATA_interp] =Determine_qinf_qsup_P1D_DD(CASES,NAMEsnap_base,DATAoffline,DATA_interp) ;





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

SNAPstressSTWOproj = cell(1,length(CASES)) ;
SNAPstressPonePROJ = cell(1,length(CASES)) ;
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
if ~isempty(BstRED_r)
    BstRED_r_normalized = BstRED_r/norm(BstRED_r,'fro'); 
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
            
            coeff =  BasisU(:,1:nREDcoor)'*(Kll*SNAP_cluster.DISP.U(DOFl,:)) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            qL  = coeff*SNAP_cluster.DISP.V' ;
            %   qL_extended = tauNON(qL) ;
            [qL_extended] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA_evaluateTAU_and_DER) ;
            qLATENT_LOC{iloc} = qL_extended(idimLAT,:) ;
        end
        
        
        
        
        % ------------------------------------------------
        %         dL = BasisU*qL_extended; % Snapshot displacements (local)
        %         dR = bsxfun(@times,SNAP_cluster.DISP.U(DOFr,:)',SNAP_cluster.DISP.S)' ;
        %         dR = dR*SNAP_cluster.DISP.V' ;
        %         ndof = size(dL,1)+size(dR,1) ;
        %         d = zeros(ndof,size(dL,2)) ;
        %         d(DOFl,:)  = dL ;
        %         d(DOFr,:)  = dR ;
        
        d = ObtainDisplacementFromSnapshotsProj(BasisU,qL_extended,SNAP_cluster,DOFr,OTHER_output,DOFl) ;
        
        
        
        %
        [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR] = SnapStressFromDispLOC...
            (VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
            SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj) ;
        
        
        %
        %
        %         % 2. Deformation gradient at all Gauss points
        %         FgradST = OPERFE.Bst*d + repmat(OPERFE.IDENTITY_F,1,size(d,2)) ;
        %         % 3. Green-Lagrante strains at all Gauss points
        %         GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
        %         % 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
        %         SNAPstressSTWOproj_LOC{iloc} = zeros(size(GLSTRAINS)) ;
        %         for isnap = 1:size(GLSTRAINS,2)
        %             [SNAPstressSTWOproj_LOC{iloc}(:,isnap) ]= PK2stress_Constitutive_Model(GLSTRAINS(:,isnap),MATPRO,DATA) ;
        %         end
        %         % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
        %         SNAPstressPonePROJ_LOC{iloc} = PK1stress(SNAPstressSTWOproj_LOC{iloc},FgradST,DATA.MESH.ndim) ;
        %
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
              BstRED_r_input = BstRED_r_normalized*norm(BstRED_l_q,'fro') ; 
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
%Ael = A_internalFORCES_ECM(:,OTHER_output.ind_elastic) ;

if DATAoffline.USE_ELEMENT_BASED_APPROACH == 1
    [wST,A_internalFORCES_ECM] = ConvertECM_points2elements(DATA,A_internalFORCES_ECM,OPERFE.wSTs)  ;
else
    wST = OPERFE.wSTs;
end


DATAoffline = DefaultField(DATAoffline,'errorFINT_linearpart',1e-6) ;

% LINEAR BASIS MATRIX
TOLloc = DATAoffline.errorFINT_linearpart;
DATAsvdLOC.HIDE_OUTPUT =  1;
[Qel_nw,SSS,VVV]  = SRSVD(A_internalFORCES_ECM(:,OTHER_output.ind_elastic),TOLloc,DATAsvdLOC)  ;
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
save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','Kll','BstRED_reactions' )




%save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','Kll')



diary off

