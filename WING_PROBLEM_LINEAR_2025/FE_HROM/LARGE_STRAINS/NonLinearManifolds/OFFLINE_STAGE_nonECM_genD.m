function OFFLINE_STAGE_nonECM_genD(DATAoffline,DATA_interp,...
    NAME_BASE,DATA_GENGAUSS,InputDataForForces)
%--------------------------------------------------------------------------
% This is a modification of 
% function OFFLINE_STAGE_nonECM_gen(DATAoffline, DATA_interp, ...
%     NAME_BASE, DATA_GENGAUSS, InputDataForForces )
%
%  whose description is shown below. The modification is related with the
%  fact that the former version was not able to deal with non-homogeneous
%  Dirichlet boundary conditions
%  Modification made on July 27th 2025, Sunday, Boraĉ, Serbia 
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/05_meta_AECM.mlx


% function OFFLINE_STAGE_nonECM_gen(DATAoffline, DATA_interp, ...
%     NAME_BASE, DATA_GENGAUSS, InputDataForForces )
% PURPOSE:
%   Executes the offline stage for constructing a **nonlinear HROM**
%   (Hyper-Reduced-Order Model) without relying on a prebuilt empirical cubature
%   mesh. It performs:
%     - Generation of reduced bases (displacements, stresses),
%     - Identification of nonlinear mapping τ(q),
%     - Derivatives τ′(q), τ″(q) (if applicable),
%     - Projection of stresses (PK2, PK1),
%     - Assembly of internal force snapshots,
%     - Point selection and weight computation via Discrete ECM (or Continuous ECM),
%     - Construction of master/slave mapping η(q) for reduced integration.
%
% INPUTS:
%   - DATAoffline         : Struct containing tolerances and flags for model reduction,
%                           stress projection, and hyperreduction control.
%   - DATA_interp         : Struct for interpolation settings of τ(q), including:
%                             • METHOD_INTERP ('BSPLINES_LEAST_SQUARES', 'NN', etc.)
%                             • METHOD_SELECT_REFERENCE_MODE ('FIRST_SNAPSHOT' or 'FIRST_SVD_MODE')
%                             • NSAMPLES, INCLUDE_SECOND_DERIVATIVES, etc.
%   - NAME_BASE           : Base name for loading training simulation results.
%   - DATA_GENGAUSS       : (Optional) Struct containing Continuous ECM options.
%   - InputDataForForces  : Function handle returning a cell array of load cases.
%
% OUTPUT (Saved in OFFLINE.mat):
%   - BasisU              : Reduced basis for displacements (restricted to DOFl).
%   - BasisStwo           : Reduced basis for PK2 stresses.
%   - ECMdata             : Struct with quadrature points, weights, and nonlinear mapping η(q).
%   - tauNON              : Nonlinear mapping τ(q) of reduced coordinates to displacement space.
%   - tauNONder           : Derivative τ′(q).
%   - tauNONder2          : Second derivative τ″(q), only if supported by method.
%   - DATA                : Updated data structure, including DOFr (constrained DOFs).
%
% KEY STEPS:
%   1. Define training cases from input loading function.
%   2. Select a reference subspace using either the first snapshot or the first SVD mode.
%   3. Project displacement snapshots into (qINF, qSUP) space and fit τ(q).
%   4. Recover PK2 stress using constitutive laws and displacement encodings.
%   5. Assemble internal force snapshots for hyperreduction.
%   6. Build a reduced-order operator for PK2 stress projection using randomized SVD.
%   7. Apply Discrete ECM (or Continuous ECM) to select optimal integration points.
%   8. Construct master/slave η mapping between ECM points and reconstruct stress integration.
%
% REMARKS:
%   - The method supports multiple interpolation strategies for τ(q), including:
%       • Shape-preserving interpolants
%       • B-splines (including least-squares fit)
%       • Shallow neural networks (NN)
%   - Second derivatives τ″(q) are included only if `INCLUDE_SECOND_DERIVATIVES = 1`.
%   - A master-slave mapping (η, η′, η″) is computed for further compression of quadrature.
%
% DEPENDENCIES:
%   - Determine_qinf_qsup.m / Determine_qinf_qsup_1SVD.m
%   - BsplinesLeastSquares.m, ShapePreInterpLOC.m, NNsINTERP.m
%   - RSVDqp.m, DiscreteECM_givenAmat.m, DiscreteECM_adaptWEIGHTS.m
%   - PK2stress_Constitutive_Model.m, PK1stress.m
%   - BasisF_from_BasisStress_PK1.m
%
% RELATED FILES:
%   - SNAPSHOTS/BEAM2D_3_param_*/INFO_SNAPSHOTS.mat
%   - MODES/*DISPmodes.msh (GiD visualization)
%   - OFFLINE.mat (output file with all HROM operators)
%
% EXAMPLE USAGE:
%   OFFLINE_STAGE_nonECM_gen(DATAoffline, DATA_interp, ...
%       'BEAM2D_3_param_', [], @InputForcesBeam)
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Comments updated by ChatGPT, 17-JUL-2025
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
    DATAoffline.errorDISP = 1e-6;
    DATAoffline.errorFINT = 1e-3;%
    DATA_interp.METHOD_SELECT_REFERENCE_MODE =  'FIRST_SVD_MODE' ; 'FIRST_SNAPSHOT';
    DATA_interp.METHOD_INTERP =  'BSPLINES_LEAST_SQUARES';  'SHAPE_PRESERVING_INTERPOLANT' ;  'NN'; 'SPLINE';    'NN' ;
    DATA_interp.NSAMPLES = 100;
    DATA_interp.order_Bslines = 4;
    
    
    
    DATA_interp.PortionExtrapolation_plot = 0 ;
    
    DATA_interp.INCLUDE_SECOND_DERIVATIVES =1;
    DATAoffline.Hyperreduction_METHOD_projection = 'STANDARD';
    DATA_interp.DATA_interp.MAKE_SVD_AMATRIX_per_project = 0 ;
    
    
    %DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints = 1;
    
    
    
    DATAoffline.USE_ELEMENT_BASED_APPROACH = 0;
    DATAoffline.errorSTRESS = 1e-2;   % For each block. Just for check that the basis matrix for displacements is representative
    DATAoffline.errorECM = 0;
    DATAoffline.errorPK2stress_basis = 1e-5;
    
    
    
    NAME_BASE = 'BEAM2D_3_param_';
    %     NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
    %     NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];
    %NAMEOFFLINEstore_CECM = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE_CECM.mat'];
    
    % Continuous ECM
    % ----------------**************************************************************************************************+
    DATA_GENGAUSS = [] ;
    
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

CASES = 1:length(InputForces) ;  % Number of training projects






switch   DATA_interp.METHOD_SELECT_REFERENCE_MODE
    case 'FIRST_SNAPSHOT'
        error('Revised this option, it has not been updated, 27th July 2025')
        [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup(CASES,NAMEsnap_base,DATAoffline) ;
    case 'FIRST_SVD_MODE'
        [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_1SVDd(CASES,NAMEsnap_base,DATAoffline) ;
end



switch DATA_interp.METHOD_INTERP
    case 'NN'
        [tauNON, tauNONder] = NNsINTERP(DATA_interp.nneurons,qINF,VV,SS,UU) ;
        %  [net] = NNsINTERP_2(DATA_interp.nneurons,qINF,VV,SS) ;
        %  error('Cannot proceed any longer .... Derivatives cannot be computed with NNs')
        nREDcoor = 1;
    case 'SPLINE'
        % DEFopt = 1 ;
        % if DEFopt == 1
        [tauNON, tauNONder,nREDcoor] = SplineInterpLOC(DATA_interp,qINF,VV',UU,SS);
    case 'SPLINE_SMOOTHING'
        % DATA_interp.NSAMPLES = 1900;
        method = 'smoothing';  % Options: 'interpolant', 'smoothing', 'least_squares'
        
        % For smoothing spline
        method_params.p = 0.95;
        
        [make_NN_evaluator, tauNONder, nREDcoor] =...
            SplineInterpLOC2(DATA_interp, qINF, VV', UU, SS, method, method_params) ;
        
        
    case 'SHAPE_PRESERVING_INTERPOLANT'
        [tauNON, tauNONder,nREDcoor] = ShapePreInterpLOC(DATA_interp,qINF,VV',UU,SS);
    case  'BSPLINES_LEAST_SQUARES'
        
        [tauNON, tauNONder,tauNONder2, nREDcoor] = BsplinesLeastSquares(DATA_interp, qINF, VV', UU, SS) ;
end

BasisU = [PhiLIN,PhiNON] ;



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
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        % Projected displacements
        % -----------------------
        if DATAoffline.UseEncoderToDetermineStresses == 0
            coeff = BasisU'*SNAP_cluster.DISP.U(DOFl,:) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            qL_extended = coeff*SNAP_cluster.DISP.V' ;
        else
            coeff =  BasisU(:,1:nREDcoor)'*SNAP_cluster.DISP.U(DOFl,:) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            qL  = coeff*SNAP_cluster.DISP.V' ;
            qL_extended = tauNON(qL) ;
        end
        
        % ------------------------------------------------
        dL = BasisU*qL_extended; % Snapshot displacements (local)
        dR = bsxfun(@times,SNAP_cluster.DISP.U(DOFr,:)',SNAP_cluster.DISP.S)' ;
        dR = dR*SNAP_cluster.DISP.V' ;
        ndof = size(dL,1)+size(dR,1) ;
        d = zeros(ndof,size(dL,2)) ;
        d(DOFl,:)  = dL ;
        d(DOFr,:)  = dR ;
        %
        % 2. Deformation gradient at all Gauss points
        FgradST = OPERFE.Bst*d + repmat(OPERFE.IDENTITY_F,1,size(d,2)) ;
        % 3. Green-Lagrante strains at all Gauss points
        GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
        % 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
        SNAPstressSTWOproj_LOC{iloc} = zeros(size(GLSTRAINS)) ;
        for isnap = 1:size(GLSTRAINS,2)
            [SNAPstressSTWOproj_LOC{iloc}(:,isnap) ]= PK2stress_Constitutive_Model(GLSTRAINS(:,isnap),MATPRO,DATA) ;
        end
        % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
        SNAPstressPonePROJ_LOC{iloc} = PK1stress(SNAPstressSTWOproj_LOC{iloc},FgradST,DATA.MESH.ndim) ;
        
        % INTERNAL FORCES MATRIX (FOR HYPERREDUCTION PURPOSES)
        A_fint = cell(1,size(qL_extended,2)) ;
        
        for itimeLOC = 1:size(qL_extended,2)
            q =  qL_extended(1:nREDcoor,itimeLOC)  ;
            tauNONder_q = tauNONder(q);
            BstRED_l_q = BstRED_l*tauNONder_q ;
            Pk1_stress = SNAPstressPonePROJ_LOC{iloc}(:,itimeLOC);
            A_fint{itimeLOC} = BasisF_from_BasisStress_PK1(BstRED_l_q,Pk1_stress,DATA);
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
    disp(['Project = ',num2str(iproj),'; ERROR stress PK2= ',num2str(STRESS_PK2_error(iproj))]);
    if STRESS_PK2_error(iproj) > DATAoffline.errorSTRESS
        %  dbstop('129')
        error('STRESS ERROR CRITERION NOT MET: TRY WITH A LOWER ERROR TOLERANCE FOR DISPLACEMENTS ')
    end
    
    
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
    [UU,SS,VV] =  RSVDT(SNAPstressSTWOproj_LOC) ;
    SNAPstressSTWOproj{iproj} = bsxfun(@times,UU',SS)' ;
    %     % PK1 stresses
    %     % ----------------
    %     SNAPstressPonePROJ_LOC = cell2mat(SNAPstressPonePROJ_LOC) ;   % exact
    %     [UU,SS,VV] =  RSVDT(SNAPstressPonePROJ_LOC) ;
    %     SNAPstressPonePROJ{iproj} = bsxfun(@times,UU',SS)' ;
    
    % internal forces for_ECM
    
    
    if DATA_interp.MAKE_SVD_AMATRIX_per_project==1
        A_internalFORCES_ECM_LOC = cell2mat(A_internalFORCES_ECM_LOC) ;
        [UU,SS,VV] =  RSVDT(A_internalFORCES_ECM_LOC) ;
        A_internalFORCES_ECM{iproj} = bsxfun(@times,UU',SS)' ;
    else
        A_internalFORCES_ECM{iproj} = horzcat(A_internalFORCES_ECM_LOC{:}) ;
    end
    
    
end

if DATA_interp.MAKE_SVD_AMATRIX_per_project == 0
    A_internalFORCES_ECM =  horzcat(A_internalFORCES_ECM{:}) ;
    
end

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

if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
    
    % ******************
    % Discrete ECM
    % ******************
    if DATA_GENGAUSS.ACTIVE == 0
        [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM_givenAmat(A_internalFORCES_ECM,DATA,OPERFE.wSTs,DATAoffline) ;
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
        
        %   DATAoffline.Hyperreduction_METHOD_projection = 'MANIFOLD';
        switch DATAoffline.Hyperreduction_METHOD_projection
            case 'MANIFOLD'
                disp('Master/Slave nonlinear mapping between ECM points ')
                % Notice that now the "master" points play the role of the
                % actual ECM points (for purposes of allocating memory for internal variables, for instance
                % What happens at slave points in terms of stresese is not of corner for the method)
                [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.etaNON,ECMdata.etaNONder,...
                    ECMdata.etaNONder2] =...
                    DiscreteECM_adaptWEIGHTS(A_internalFORCES_ECM,setPoints,wRED,DATA_interp,OPERFE.wSTs) ;
                ECMdata.setElements = large2smallREP(ECMdata.setPoints,DATA.MESH.ngaus) ;
                disp(['Master element(s) =',num2str(ECMdata.setElements(:)')])
                ECMdata.setElements_slv = large2smallREP(ECMdata.setPoints_slv,DATA.MESH.ngaus) ;
                disp(['Slave element(s) =',num2str(ECMdata.setElements_slv(:)')])
            otherwise
                ECMdata.etaNON =[] ;
        end
        
        
        
    else
        error('Option not available (2-July-2025)')
        % ****************************************************************
        % Function for computing the position of the integration points
        % ****************************************************************
        Nst = OTHER_output.Nst ;
        NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
        if  ~exist(NAME_MODES_FOLDER)
            mkdir(NAME_MODES_FOLDER)
        end
        DATA_GENGAUSS.NameFileMesh_ECM = [NAME_MODES_FOLDER,NAME_BASE,'DECMpoints'] ;
        DATALOC = [] ;
        
        DATA_GENGAUSS.NameFileMesh_FINT = [NAME_MODES_FOLDER,NAME_BASE,'InternalForceModes'] ;
        DATALOC = [] ;
        
        DATA_GENGAUSS.NameFileMesh_CECM = [NAME_MODES_FOLDER,NAME_BASE,'CECMpoints'] ;
        DATALOC = [] ;
        
        
        
        [ECMdata] = ContinuousECM(BstRED_l,BasisPone,DATA,OPERFE.wSTs,DATAoffline,DATA_GENGAUSS,...
            MESH,Nst) ;
    end
    
else
    error('Option not maintained (2-Jul-2025)')
    SNAPredFINT = BasisF_from_BasisStress_PK1_ELEMS(BstRED_l,BasisPone,DATA, OPERFE.wSTs)  ;
    wSTs_LOC = ones(size(SNAPredFINT,1),1) ;
    %     sqrt_wST = sqrt(OPERFE.wSTs) ;
    %     SNAPredFINT = bsxfun(@times,SNAPredFINT,sqrt_wST) ;
    % Determine an  orthogonal basis matrix $Q$ for the column space of $\SNAPredFINTw{}{}$
    %DATAoffline.errorFINT = 1e-3;
    DATAsvd.RELATIVE_SVD = 1;
    [Q,S,V,eSVD,Rsup] = RSVDT(SNAPredFINT,DATAoffline.errorFINT,[],0,DATAsvd) ;
    
    if DATAoffline.errorFINT == 0
        ifig = 3000 ;
        SVDplotERROR_local(S,ifig) ;
    end
    
    % % Enlarge the basis matris for SNAPredFINT
    a  = wSTs_LOC - Q*(Q'*wSTs_LOC) ;
    if norm(a) > 1e-10
        a = a/norm(a) ;
        Q = [Q,a] ;
    end
    % Empirical cubature method
    % -------------------------
    DATA_ECM = [] ;
    DATA_ECM.TOL = DATAoffline.errorECM ;
    [setElements,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs_LOC,DATA_ECM)  ;
    disp(['Element-based approach *************************++'])
    disp(['List of selected m = ',num2str(length(setElements)),' elements'])
    disp(num2str(setElements'))
    figure(345)
    hold on
    xlabel('Points')
    ylabel('Weights')
    
    bar(sort(wRED,'descend'))
    
    % Determine set of points
    setPoints = small2large(setElements,DATA.MESH.ngaus_STRESS) ;
    disp(['Total number of Gauss points = ',num2str(length(setPoints))])
    wRED = repmat(wRED',DATA.MESH.ngaus_STRESS,1) ;
    wRED = wRED(:).*OPERFE.wSTs(setPoints,:) ;
    
    ECMdata.setPoints = setPoints ;
    ECMdata.wRED = wRED ;
    ECMdata.setElements = setElements ;
    
end

disp(['*********************************************************************'])
disp(['Number of displacement modes = ',num2str(size(BasisU,2))])
disp(['Number of PK2-stress modes = ',num2str(size(BasisStwo,2))])
%disp(['Number of PK1-stress modes = ',num2str(size(BasisPone,2))])



% STORING INFORMATION ---HYPERREDUCED-ORDER OPERATORS

BASES.BasisU = BasisU ;
DATA.DOFr = DOFr;
save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','tauNON','tauNONder','tauNONder2')



diary off

