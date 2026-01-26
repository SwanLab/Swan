function OFFLINE_STAGE_nonECM_plastK_ecmADAPT(DATAoffline,DATA_interp,...
    NAME_BASE,DATA_GENGAUSS,InputDataForForces)
% Modification of OFFLINE_STAGE_nonECM_plastK, for adaptive ECM without
% separation of linear and nonlinear components
% 9-SEpt-2025, Balmes 185, BArcelona
% =========================================================================
% OFFLINE STAGE — NONLINEAR MANIFOLD HROM (PLASTICITY) WITH Kll METRIC
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
if isempty(DATA_GENGAUSS)
    % Continuous ECM
    % ----------------**************************************************************************************************+
    DATA_GENGAUSS.ACTIVE =0;
    DATA_GENGAUSS.ENFORCE_SUM_WEIGHTS_EQUAL_VOLUME =1;
    
    DATA_GENGAUSS.TOL =  1e-8 ; % 1e-10 ; % Tolerance removing points
    DATA_GENGAUSS.TOL_low = 1e-6 ; % Tolerance removing points
    DATA_GENGAUSS.kMAX = 40 ;
    DATA_GENGAUSS.PLOT_INTERNAL_FORCE_MODES = 1;
    DATA_GENGAUSS.PLOT_REDUCED_GAUSS_IN_GID =    0;  % Print location POINTS in GID MESH
    DATA_GENGAUSS.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS = 1  ; % Factor printing points (0 <TOL <1)
    DATA_GENGAUSS.RESTRICTED_DOMAIN_SELECTION_BYLAYERS =  1;  % Number of layers to be excluded for ECM points
    DATA_GENGAUSS.WHICH_BOUNDARIES_TO_EXCLUDE =  [ ] ; [1,3:4]; % 3:6;  % bOUNDARIES OF ELEMENTS ARE TO BE EXCLUDED.
end



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
    DATA_evaluateTAU_and_DER, nREDcoor,Kll] =Determine_qinf_qsup_PLAST_1DK(CASES,NAMEsnap_base,DATAoffline,DATA_interp) ;





BasisU = [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl] ;



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

SNAPstressSTWOproj = cell(1,length(CASES)) ;
SNAPstressPonePROJ = cell(1,length(CASES)) ;
A_internalFORCES_ECM = cell(1,length(CASES)) ; % For hyperreduction purposes
BstRED_l = OPERFE.Bst(:,DOFl)*BasisU ;
qPLAST_snap =  cell(1,length(CASES));

DATAoffline = DefaultField(DATAoffline,'UseEncoderToDetermineStresses',0);
%DATAoffline.UseEncoderToDetermineStresses = 1;

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
      qPLAST_loc = cell(1,length(NAME_SNAP_loc)); % A_fint per loc
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        % Projected displacements
        % -----------------------
%         DATAoffline.UseEncoderToDetermineStresses = 1;
%         if DATAoffline.UseEncoderToDetermineStresses == 0
%             
%             coeff = BasisU'*(Kll*SNAP_cluster.DISP.U(DOFl,:)) ;
%             coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
%             qL_extended = coeff*SNAP_cluster.DISP.V' ;
%         else
            
            coeff =  BasisU(:,1:nREDcoor)'*(Kll*SNAP_cluster.DISP.U(DOFl,:)) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            qL  = coeff*SNAP_cluster.DISP.V' ;
            %   qL_extended = tauNON(qL) ;
            [qL_extended] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA_evaluateTAU_and_DER) ;
            ind_latent_plastic = 2;
        qPLAST_loc{iloc} = qL(ind_latent_plastic,:) ;
      %  end
        
        % ------------------------------------------------
        dL = BasisU*qL_extended; % Snapshot displacements (local)
        dR = bsxfun(@times,SNAP_cluster.DISP.U(DOFr,:)',SNAP_cluster.DISP.S)' ;
        dR = dR*SNAP_cluster.DISP.V' ;
        ndof = size(dL,1)+size(dR,1) ;
        d = zeros(ndof,size(dL,2)) ;
        d(DOFl,:)  = dL ;
        d(DOFr,:)  = dR ;
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
            A_fint{itimeLOC} = BasisF_from_BasisStress_PK1(BstRED_l_q,Pk1_stress,DATA);
        end
        
        % DATA_interp = DefaultField(DATA_interp,'MAKE_SVD_AMATRIX_per_project',1) ;
        
        %  if DATA_interp.MAKE_SVD_AMATRIX_per_project==1
        
        %      A_internalFORCES_ECM_LOC{iloc} = cell2mat(A_fint) ;
        % else
        A_internalFORCES_ECM_LOC{iloc} = A_fint ;
        % end
        
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
    RELTOL = 1e-4;  % WE only need for reconstruction
    [UU,SS,VV] =  SRSVD(SNAPstressSTWOproj_LOC,RELTOL) ;
    SNAPstressSTWOproj{iproj} = bsxfun(@times,UU',SS)' ;
    %     % PK1 stresses
    %     % ----------------
    %     SNAPstressPonePROJ_LOC = cell2mat(SNAPstressPonePROJ_LOC) ;   % exact
    %     [UU,SS,VV] =  RSVDT(SNAPstressPonePROJ_LOC) ;
    %     SNAPstressPonePROJ{iproj} = bsxfun(@times,UU',SS)' ;
    
    % internal forces for_ECM
    
    
    %     if DATA_interp.MAKE_SVD_AMATRIX_per_project==1
    %         A_internalFORCES_ECM_LOC = cell2mat(A_internalFORCES_ECM_LOC) ;
    %       %  [UU,SS,VV] =  RSVDT(A_internalFORCES_ECM_LOC) ;
    %         A_internalFORCES_ECM{iproj} = A_internalFORCES_ECM_LOC ;
    %         % compress the information of the internal forces
    %         % Note that what counts is the column space of this matrix
    %     else
    A_internalFORCES_ECM{iproj} = horzcat(A_internalFORCES_ECM_LOC{:}) ;
      qPLAST_snap{iproj} = cell2mat(qPLAST_loc);
    % end
    
    
end

%if DATA_interp.MAKE_SVD_AMATRIX_per_project == 0
A_internalFORCES_ECM =  horzcat(A_internalFORCES_ECM{:}) ;
qPLAST_snap = cell2mat(qPLAST_snap) ;
 qPLAST= qPLAST_snap(:,OTHER_output.ind_plastic)  ;
 % In the linear range, snaphots do not change from cluster to cluster, so
 % we make
 
 A_internalFORCES_ECM_lin = (A_internalFORCES_ECM(:,OTHER_output.ind_elastic(end))) ; 
 % So now we make the following: 
 A_internalFORCES_ECM_non = cell(1,length(OTHER_output.ind_plastic)+1) ; 
 A_internalFORCES_ECM_non{1} = A_internalFORCES_ECM_lin{1} ; 
 A_internalFORCES_ECM_non(2:end) = A_internalFORCES_ECM(:,OTHER_output.ind_plastic) ; 
 qPLAST = [0,qPLAST(:)'] ;

% BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK2 STRESSES (FOR RECONSTRUCTION PURPOSES)
% -------------------------------------------------------------------------------------------------
TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAPstressSTWOproj),1)' ;

DATAsvd=[];
[BasisStwo,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressSTWOproj,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
disp('***********************************************************')
%
% 
% -------------------
% HYPERREDUCTION
% -------------------
  ECMdata = SAW_ECM_elastplastJOIN(DATAoffline,A_internalFORCES_ECM_non,DATA,OPERFE,qPLAST)  ;
% 
% % ******************
% % Discrete ECM
% % ******************
% %if DATA_GENGAUSS.ACTIVE == 0
% [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM_givenAmat(A_internalFORCES_ECM,DATA,OPERFE.wSTs,DATAoffline) ;
% ECMdata.setPoints = setPoints ;
% ECMdata.wRED = wRED ;
% proporP = length(setPoints)/DATA.MESH.ngausT*100;
% disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])
% 
% setElements = large2smallREP(setPoints,DATA.MESH.ngaus) ;
% disp('****************************+')
% disp(['List of selected m = ',num2str(length(setElements)),' elements'])
% disp(num2str(setElements'))
% clipboard('copy',num2str(setElements'));
% ECMdata.setElements = setElements ;
% 
% %   DATAoffline.Hyperreduction_METHOD_projection = 'MANIFOLD';
% switch DATAoffline.Hyperreduction_METHOD_projection
%     case 'MANIFOLD'
%         disp('Master/Slave nonlinear mapping between ECM points ')
%         % Notice that now the "master" points play the role of the
%         % actual ECM points (for purposes of allocating memory for internal variables, for instance
%         % What happens at slave points in terms of stresese is not of corner for the method)
%         
%         DATA_interp =DefaultField(DATA_interp,'NumberOfMasterPOINTS_ECMadaptive',1);
%         
%         if DATA_interp.NumberOfMasterPOINTS_ECMadaptive == 1
%             [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.etaNON,ECMdata.etaNONder,...
%                 ECMdata.etaNONder2] =...
%                 DiscreteECM_adaptWEIGHTS(A_internalFORCES_ECM,setPoints,wRED,DATA_interp,OPERFE.wSTs) ;
%         elseif DATA_interp.NumberOfMasterPOINTS_ECMadaptive == 2
%             [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.etaNON,ECMdata.etaNONder,...
%                 ECMdata.etaNONder2] =...
%                 DiscreteECM_adaptWEIGHTS_2p(A_internalFORCES_ECM,setPoints,wRED,DATA_interp,OPERFE.wSTs) ;
%         else
%             error('General case not implemented yet (5th August 2025)')
%         end
%         
%         
%         
%         
%         ECMdata.setElements = large2smallREP(ECMdata.setPoints,DATA.MESH.ngaus) ;
%         disp(['Master element(s) =',num2str(ECMdata.setElements(:)')])
%         ECMdata.setElements_slv = large2smallREP(ECMdata.setPoints_slv,DATA.MESH.ngaus) ;
%         disp(['Slave element(s) =',num2str(ECMdata.setElements_slv(:)')])
%     otherwise
%         ECMdata.etaNON =[] ;
% end


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
%
% end

disp(['*********************************************************************'])
disp(['Number of displacement modes = ',num2str(size(BasisU,2))])
disp(['Number of PK2-stress modes = ',num2str(size(BasisStwo,2))])
%disp(['Number of PK1-stress modes = ',num2str(size(BasisPone,2))])



% STORING INFORMATION ---HYPERREDUCED-ORDER OPERATORS

%BASES.BasisU = BasisU ;
DATA.DOFr = DOFr;
DATA.DATA_evaluateTAU_and_DER = DATA_evaluateTAU_and_DER;
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
save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','Kll' )




%save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','Kll')



diary off

