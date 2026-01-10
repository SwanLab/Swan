function OFFLINE_STAGE_nonECM_genFASTd(DATAoffline,DATA_interp,...
    NAME_BASE,DATA_GENGAUSS,InputDataForForces,DATA_interp_ECM)
%--------------------------------------------------------------------------
% OFFLINE_STAGE_nonECM_genFASTd
%
% PURPOSE
%   Build the OFFLINE dataset for a nonlinear HROM **without** a predefined
%   empirical–cubature mesh, adapted to **non-homogeneous Dirichlet BCs** and
%   optionally including **boundary-reaction modes** in the force snapshot set.
%   The routine trains an encoder τ(q) (and its derivatives) on displacement
%   snapshots, reconstructs stresses, assembles internal-force columns, and
%   performs (Discrete/SAW) ECM selection. It saves everything needed for the
%   ONLINE stage into OFFLINE.mat.
%
% MATHEMATICAL/MECHANICAL CONTEXT (tie-in with the slides)
%   • Weak form & virtual work. Internal forces are assembled from the weak form
%     (virtual work) and, for nonlinear problems, solved via Newton–Raphson with
%     tangent K = ∂F̂/∂d. These ideas motivate how we form force snapshots and
%     their Jacobian surrogates offline. :contentReference[oaicite:0]{index=0}
%   • Vectorized global operators. We exploit precomputed sparse operators to
%     assemble internal forces in matrix form:
%         F̂ = Bᵀ ω S  ,   E = B d  ,   U = N d ,
%     where L, B, N and the quadrature weights ω come from the FE assembly on
%     Gauss points (globalized Boolean/shape-function operators). This matches
%     the vectorized formulation shown in the slides (Eqs. (160)–(174)). :contentReference[oaicite:1]{index=1}
%   • Gauss quadrature. Element contributions to forces/stiffness use standard
%     Gauss rules (p ≤ 2m−1 exactness) and the parent→physical mapping with
%     he/2 scaling; we respect that when building reduced force columns. :contentReference[oaicite:2]{index=2}
%
% MANIFOLD HROM VIEW (what “manifold” adds)
%   • Displacements are represented as d ≈ Φ_lin q_lin ⊕ Φ_non τ(q), where τ(q)
%     is learned (B-spline LS here) and evaluated with its derivatives τ′(q)
%     (and optionally τ″). The encoder provides consistent reduced kinematics for
%     stress/force reconstruction on the manifold.
%   • We assemble columns of the **internal-force snapshot operator** using the
%     chain rule with τ′(q): for each snapshot state,
%         Bst_red(q) = Bst * (A*BasisU) * τ′(q)      (affine BCs via A if present),
%         A_fint(:,t) = BasisF_from_BasisStress_PK1( Bst_red(q_t) , P(q_t) , DATA ).
%     If enabled, we **augment** with boundary-reaction modes so hyper-reduction
%     “sees” reaction-driven work as well (see KW:ReactF section below).
%
% NON-HOMOGENEOUS DIRICHLET BCs (affine handling)
%   • The DOFs are split as free DOFl and prescribed DOFr. We reconstruct
%     total displacements d = [dL; dR] so that F = I + Bst*d is consistent at
%     Gauss points when computing PK2/PK1 stresses (weak-form consistent). This
%     follows the slides’ treatment of essential vs. natural BCs in weak form. :contentReference[oaicite:3]{index=3}
%
% INPUTS
%   DATAoffline
%     .errorDISP, .errorFINT, .errorSTRESS, .errorPK2stress_basis
%     .Hyperreduction_METHOD_projection : 'STANDARD' | 'MANIFOLD' | 'CECM_BASED_STRATEGY' | 'SAW_ECM'
%     .USE_ELEMENT_BASED_APPROACH       : (kept = 0)
%     .UseEncoderToDetermineStresses    : use τ(q) to reconstruct dL (default 1)
%     .ECMforReactionForces             : true → include reaction modes (see below)
%     .Hyperreduction_Separate_Slave_contribution : build master/slave maps for MANIFOLD
%     .Version_MAW_ECM_withActiveSET    : plasticity-oriented CECM variant
%   DATA_interp
%     .METHOD_SELECT_REFERENCE_MODE, .METHOD_INTERP='BSPLINES_LEAST_SQUARES'
%     .NSAMPLES, .order_Bslines, .INCLUDE_SECOND_DERIVATIVES, .MAKE_SVD_AMATRIX_per_project
%     .IND_PROJECT_LINEAR (optional indices treated as linear reference only)
%   NAME_BASE
%     Base name for ./SNAPSHOTS/<NAME_BASE*> folders.
%   DATA_GENGAUSS
%     (Unused here; continuous ECM path kept for compatibility.)
%   InputDataForForces
%     Handle returning load definitions for all training projects.
%   DATA_interp_ECM
%     Settings for manifold η(q) regression when building master/slave maps.
%
% OUTPUTS / SIDE EFFECTS
%   • Saves ./SNAPSHOTS/<NAME_BASE>/OFFLINE.mat containing:
%       BasisU, BasisStwo, ECMdata, DATA, DATAoffline, BstRED_reactions
%   • Exports GiD mode files in ./MODES/.
%   • Writes a short log to OFFLINE.txt.
%
% WORKFLOW (high-level)
%   (1) Load snapshot meta-data per project; build Φ_lin (and Φ_non if used).
%   (2) Train B-spline LS encoder; get evaluators for τ, τ′ (τ″ if enabled).
%   (3) For all snapshots, reconstruct d = [dL; dR], compute F, E_GL, PK2, PK1.
%   (4) Validate PK2 reconstruction vs. exact (rel. Frobenius ≤ errorSTRESS).
%   (5) Assemble internal-force columns using Bst_red(q_t) and PK1 stresses.
%       If ECMforReactionForces==true, **augment** columns with reaction modes.
%   (6) Stack A_fint (or SVD-compress per project), then run ECM/SAW-ECM.
%   (7) If MANIFOLD hyper-reduction is selected, build η(q) master/slave maps and
%       store η′ regression data for ONLINE.
%   (8) Build PK2 stress basis (RSVDqp) for stress reconstruction ONLINE.
%
% KW:ReactF — INCLUDING BOUNDARY REACTIONS IN HYPER-REDUCTION
%   • DetermineModeReactionsOUTPUT_glo() returns a global reaction operator
%     BstRED_r. When DATAoffline.ECMforReactionForces == true we assemble
%     columns as:
%         A_fint(:,t) = BasisF_from_BasisStress_PK1( [BstRED_l*τ′(q_t), BstRED_r], P(q_t), DATA )
%     so the ECM/SAW-ECM selection accounts for volumetric + reaction work.
%
% NEWTON–RAPHSON CONSISTENCY (ties to slides)
%   • The slides derive K = ∂F̂/∂d and its vectorized form K = Bᵀ ω C B for
%     standard materials. Our offline columns mimic the same linearization path
%     (chain rule via τ′(q)) so that ONLINE tangent surrogates remain consistent
%     with the weak form and Gauss integration (cf. the derivations of the weak
%     form, Newton–Raphson, and vectorization). :contentReference[oaicite:4]{index=4}
%
% KEY FILES / DEPENDENCIES
%   Determine_qinf_qsup_LINd, Determine_qinf_qsup_1SVDd
%   BsplinesLeastSquares_fast
%   StrainGreenLagrange, PK2stress_Constitutive_Model, PK1stress
%   BasisF_from_BasisStress_PK1
%   RSVDT, RSVDqp
%   DiscreteECM_givenAmat, DiscreteECM_adaptWEIGHTSfst / fstNS, SAW_ECM_large1param
%   DetermineModeReactionsOUTPUT_glo
%
% IMPORTANT FLAGS & EFFECTS
%   • DATAoffline.ECMforReactionForces = true
%       Adds BstRED_r to force columns → ECM respects reaction-driven work.
%   • DATA_interp.MAKE_SVD_AMATRIX_per_project = 1
%       SVD-compress A_fint per project to reduce memory before concatenation.
%
% FAILURE MODES / TIPS
%   • Stress reconstruction too inaccurate → enrich BasisU (loosen errorDISP),
%     refine τ(q) (more samples / higher order / include τ″), or increase PK2 basis.
%   • Ill-conditioning in τ′(q) → regularize B-spline LS or reduce latent dim.
%
% VERSION HISTORY
%   • 2025-09-11  Added non-homog. Dirichlet handling (affine BCs).
%   • 2025-10-02  (KW:ReactF) Optional inclusion of boundary-reaction modes
%                 in A_fint and storage of BstRED_reactions.
%
% AUTHORSHIP
%   Original code: Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Comments updated: ChatGPT (2025-10-08)
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
 
InputForces = feval(InputDataForForces) ;

CASES = 1:length(InputForces.INPUTS_PARAMETERS) ;  % Number of training projects
DATA_interp = DefaultField(DATA_interp,'IND_PROJECT_LINEAR',[]);
%Determine_DATA_evaluateTAU_and_DER = true; 
DATA_interp.IgnoreFirstSTep = 0 ; 
switch   DATA_interp.METHOD_SELECT_REFERENCE_MODE
    case 'MACRO_STRAINS'
        error('Option not fully developed...')
    case 'LINEAR_MODES'
        [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_LINd(CASES,NAMEsnap_base,DATAoffline,DATA_interp) ;
    case 'FIRST_SVD_MODE'
        [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_1SVDd(CASES,NAMEsnap_base,DATAoffline) ;
        
       case 'TRANSLATIONAL_BUMP'
            warning('Tentative, not finished')
        [ PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_TRANSLATION(CASES,NAMEsnap_base,DATAoffline,NAME_BASE) ;    
        
      case 'MAXIMUM_DISPLACEMENT_LOCATION'
             
        [ PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output,DATAoffline]...
            = Determine_qinf_qsup_LocMaxDisp(CASES,NAMEsnap_base,DATAoffline,NAME_BASE) ;       
        case 'DISPLACEMENT_ONE_NODE'
             
        [ PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output,DATAoffline]...
            = Determine_qinf_qSUP_1NODE(CASES,NAMEsnap_base,DATAoffline,NAME_BASE,DATA_interp) ;    
   
     case 'FIRST_2SVD_modes'
         warning('Tentative, not finished')
        [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_2SVDd(CASES,NAMEsnap_base,DATAoffline,NAME_BASE) ;    
        error('')
        case 'ConcoctedMode'
         warning('Tentative, not finished')
        [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_CONCOCTED(CASES,NAMEsnap_base,DATAoffline,NAME_BASE) ;    
        
    case 'SEARCH_MODE_WITHIN_SPAN_RIGHT_SINGULAR_VECTORS'
        % Tempative 
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/18_movingLOAD.mlx
           [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_1SVD_search(CASES,NAMEsnap_base,DATAoffline) ;
        
        
        
    case 'STANDARD_ROM'
        qSUP = [] ; qINF=[] ; 
        [PhiLIN,PhiNON,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
            = Determine_qinf_qsup_STDROM(CASES,NAMEsnap_base,DATAoffline) ;
        
    otherwise
        error('Option not implemented')
end

BasisU = [PhiLIN,PhiNON] ;


if ~isempty(qSUP) 
    
    DATA_interp = DefaultField(DATA_interp,'SubSampling_qINELASTIC_master',false) ;
    
    if ~DATA_interp.SubSampling_qINELASTIC_master
        
        switch  DATA_interp.METHOD_SELECT_REFERENCE_MODE
            case  {'MAXIMUM_DISPLACEMENT_LOCATION','DISPLACEMENT_ONE_NODE'}
        [DATA_evaluateTAU_and_DER,nREDcoor]= BsplinesLeastSquares_MAXLOC(DATA_interp, qINF, VV', UU, SS) ; 
%             case 'DISPLACEMENT_ONE_NODE'
%                  [DATA_evaluateTAU_and_DER,nREDcoor]= BsplinesLeastSquares_DISP1node(DATA_interp, qINF, VV', UU, SS) ; 
            otherwise
            [DATA_evaluateTAU_and_DER,nREDcoor]= BsplinesLeastSquares_fast(DATA_interp, qINF, VV', UU, SS) ;     
        end
        
    else
        [DATA_evaluateTAU_and_DER,DATA_interp,nREDcoor]= BsplinesSUBSAMPL_large(qINF,qSUP,DATA_interp ) ;
     end
    
else
    %[DATA_evaluateTAU_and_DER,nREDcoor]= StandardROM(DATA_interp,BasisU ) ;
    
    nREDcoor = size(BasisU,2); 
    DATA_evaluateTAU_and_DER.nameFunctionEvaluate = 'tauFUN_identity';
    

end
%end
DATA.nREDcoor = nREDcoor ; 

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
% Now we have to compute the PK2 stresses produced by these strains, and check
% that all blocks are approximated to an accuracy theshold
% DATAoffline.\errorSTRESS
DATA_interp = DefaultField(DATA_interp,'IND_PROJECT_LINEAR',[]);

if ~isempty(DATA_interp.IND_PROJECT_LINEAR)
    CASES_orth = setdiff(CASES,DATA_interp.IND_PROJECT_LINEAR) ;
else
    CASES_orth = CASES;
end

STRESS_PK2_error = zeros(length(CASES_orth),1) ;

SNAPstressSTWOproj = cell(1,length(CASES_orth)) ;

OTHER_output = DefaultField(OTHER_output,'DISP_CONDITIONS',[]) ;
OTHER_output.DISP_CONDITIONS = DefaultField(OTHER_output.DISP_CONDITIONS,'A',[]) ;

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
% -------------
% ------------------------------------------------------------------------------------------



IndexMatrixECM.fint = 1:nREDcoor; 
if ~isempty(BstRED_r)
IndexMatrixECM.react = IndexMatrixECM.fint+1:IndexMatrixECM.fint+size(BstRED_r,2) ; 
end




%DATAoffline = DefaultField(DATAoffline,'UseEncoderToDetermineStresses',1);
DATAoffline.UseEncoderToDetermineStresses = 1;
qLATENT = cell(1,length(CASES_orth)) ;
idimLAT = 1;

DATAoffline = DefaultField(DATAoffline,'ECMforReactionForces',false) ; %   = true ;
DATAoffline = DefaultField(DATAoffline,'errorFINT_reactions',[]) ;
 A_internalFORCES_ECM = cell(1,length(CASES_orth)) ; % For hyperreduction purposes


for iproj = 1:length(CASES_orth)
    iprojLOC = CASES_orth(iproj)  ;
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iprojLOC))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    SNAPstressSTWOproj_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SNAPstressSTWO_LOC  = cell(1,length(NAME_SNAP_loc)) ;
    
    SNAPstressPonePROJ_LOC = cell(1,length(NAME_SNAP_loc)) ;
    A_internalFORCES_ECM_LOC = cell(1,length(NAME_SNAP_loc)) ;
 
    qLATENT_LOC = cell(1,length(NAME_SNAP_loc)) ;
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        % Projected displacements
        % -----------------------
        if DATAoffline.UseEncoderToDetermineStresses == 0
            % Just for verification purposes, not used on a regular basis
            coeff = BasisU'*SNAP_cluster.DISP.U(DOFl,:) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            qL_extended = coeff*SNAP_cluster.DISP.V' ;
            tauNONder = [] ;
        else
            
            switch DATA_interp.METHOD_SELECT_REFERENCE_MODE
            case  'MAXIMUM_DISPLACEMENT_LOCATION'
                [qLATENT_LOC,qL_extended] = EncoderMaxDISP(DATAoffline,SNAP_cluster,DOFl,DATA_evaluateTAU_and_DER,qLATENT_LOC,iloc,MESH) ; 
               tauNONder = [] ; 
                case 'DISPLACEMENT_ONE_NODE'
                    [qLATENT_LOC,qL_extended] = Encoder1Node(DATAoffline,SNAP_cluster,DOFl,DATA_evaluateTAU_and_DER,qLATENT_LOC,iloc,MESH,...
                        DATA_interp) ; 
               tauNONder = [] ; 
                otherwise
            [qLATENT_LOC,qL_extended] = Encoder1MasterMode(BasisU,nREDcoor,SNAP_cluster,DOFl,DATA_evaluateTAU_and_DER,qLATENT_LOC,iloc,idimLAT) ; 
          tauNONder = [] ; 
            end
            
            
            
        end
        
        % ------------------------------------------------
        
        d = ObtainDisplacementFromSnapshotsProj(BasisU,qL_extended,SNAP_cluster,DOFr,OTHER_output,DOFl) ;        
        
        %
        % 2. Deformation gradient at all Gauss points
        FgradST = OPERFE.Bst*d + repmat(OPERFE.IDENTITY_F,1,size(d,2)) ;
        % 3. Green-Lagrante strains at all Gauss points
        GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
        % 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
        SNAPstressSTWOproj_LOC{iloc} = zeros(size(GLSTRAINS)) ;
        for isnap = 1:size(GLSTRAINS,2)
            [SNAPstressSTWOproj_LOC{iloc}(:,isnap) ]= PK2stress_Constitutive_Model(GLSTRAINS(:,isnap),MATPRO,DATA,FgradST(:,isnap)) ;
        end
        % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
        SNAPstressPonePROJ_LOC{iloc} = PK1stress(SNAPstressSTWOproj_LOC{iloc},FgradST,DATA.MESH.ndim) ;
        
        % INTERNAL FORCES MATRIX (FOR HYPERREDUCTION PURPOSES)
        A_fint = cell(1,size(qL_extended,2)) ;
         for itimeLOC = 1:size(qL_extended,2)
            q =  qL_extended(1:nREDcoor,itimeLOC)  ;
            if isempty(tauNONder)
                [~,tauNONder_q,~] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,q,DATA_evaluateTAU_and_DER) ;
            else
                tauNONder_q = tauNONder(:,itimeLOC) ;
            end
            %    tauNONder_q = tauNONder(q);
            
            BstRED_l_q = BstRED_l*tauNONder_q ;
            Pk1_stress = SNAPstressPonePROJ_LOC{iloc}(:,itimeLOC);
            
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
          %  A_react_ECM_LOC{iloc} = cell2mat(A_fint) ;
        else
            A_internalFORCES_ECM_LOC{iloc} = A_fint ;
          %  A_react_ECM_LOC{iloc} = A_react ;
        end
        
    end
    
    
    SNAPstressSTWOproj_LOC = cell2mat(SNAPstressSTWOproj_LOC) ;   % Approximate
    SNAPstressSTWO_LOC = cell2mat(SNAPstressSTWO_LOC) ;   % exact
    
    % Check if it meets the error criterion
    if  size(SNAPstressSTWO_LOC,2) == size(SNAPstressSTWOproj_LOC,2)
    STRESS_PK2_error(iproj) = norm(SNAPstressSTWO_LOC-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;
    else
        STRESS_PK2_error(iproj) = norm(SNAPstressSTWO_LOC(:,2:end)-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;  
    end
    disp(['Project = ',num2str(iproj),'; ERROR stress PK2= ',num2str(STRESS_PK2_error(iproj))]);
    if STRESS_PK2_error(iproj) > DATAoffline.errorSTRESS
        %  dbstop('129')
        error('STRESS ERROR CRITERION NOT MET: TRY WITH A LOWER ERROR TOLERANCE FOR DISPLACEMENTS ')
    end
    
    
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
   % [UU,SS,VV] =  RSVDT(SNAPstressSTWOproj_LOC) ;
    SNAPstressSTWOproj{iproj} = [] ; %bsxfun(@times,UU',SS)' ;
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
      %  A_react_ECM{iproj} = horzcat(A_react_ECM_LOC{:}) ;
    end
    qLATENT{iproj} = cell2mat(qLATENT_LOC) ;
    
end

if DATA_interp.MAKE_SVD_AMATRIX_per_project == 0
    A_internalFORCES_ECM =  horzcat(A_internalFORCES_ECM{:}) ;
   %  A_react_ECM =  horzcat(A_react_ECM{:}) ;
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
% if DATAoffline.USE_ELEMENT_BASED_APPROACH == 0
%
%     % ******************
%     % Discrete ECM
%     % ******************
%     if DATA_GENGAUSS.ACTIVE == 0

switch DATAoffline.Hyperreduction_METHOD_projection
    
    case 'ALL_GAUSS_POINTS'
     ECMdata.setPoints = 1:length(OPERFE.wSTs) ;
        ECMdata.wRED = OPERFE.wSTs ;
         disp('All Gauss points selected')
        
        setElements = 1:DATA.MESH.nelem;  ;
        
       % clipboard('copy',num2str(setElements'));
        ECMdata.setElements = setElements ;
    case  'SAW_ECM'
        error('Option no longer used....')
        %ECMdata = SAW_ECM_large1param(DATAoffline,A_internalFORCES_ECM,DATA,OPERFE,qLATENT,DATA_interp)  ;
    otherwise
        DATA.IndexMatrixECM = IndexMatrixECM; 
        [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM_givenAmat(A_internalFORCES_ECM,DATA,OPERFE.wSTs,DATAoffline) ;
        ECMdata.setPoints = setPoints ;
        ECMdata.wRED = wRED ;
        proporP = length(setPoints)/DATA.MESH.ngausT*100;
        disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])
        
        setElements = large2smallREP(setPoints,DATA.MESH.ngaus) ;
        disp('****************************+')
        disp(['List of selected m = ',num2str(length(setElements)),' elements'])
        disp(num2str(setElements'))
       % clipboard('copy',num2str(setElements'));
        ECMdata.setElements = setElements ;
        
        %   DATAoffline.Hyperreduction_METHOD_projection = 'MANIFOLD';
        switch DATAoffline.Hyperreduction_METHOD_projection
            case 'CECM_BASED_STRATEGY'
                disp('Hyperreduction strategy based on the continuous ECM ')
                
                
                DATAoffline = DefaultField(DATAoffline,'Version_MAW_ECM_withActiveSET',false) ;
                
                DATAoffline.Index_q_used_for_regression_WEIGHTS_SAW_ECM = 1;
                
                if DATAoffline.Version_MAW_ECM_withActiveSET
                    % Version devised for plasticity problems, adapted on
                    % 30-Sept-2025
                    % SEe /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
                    LOCALdata = [] ; Uel = [] ;
                    ECMdata = CECM_based_ManifAdWeights_PLaCE(A_internalFORCES_ECM,DATA,OPERFE.wSTs,DATA_interp,setPoints,wRED,...
                        qLATENT,LOCALdata,Uel,DATAOUT.Qbasis_weighted,DATAoffline) ;
                    
                else
                    % Version before 30-Sept-2025
                    %   devised on 17th September 2025, see
                    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
                    ECMdata = CECM_based_ManifAdWeights(A_internalFORCES_ECM,DATA,OPERFE.wSTs,DATA_interp,setPoints,wRED,qLATENT) ;
                    
                end
                
                
                
            case 'MANIFOLD'
                error('Option no longer used...')
                %                 disp('Master/Slave nonlinear mapping between ECM points ')
                %                 % Notice that now the "master" points play the role of the
                %                 % actual ECM points (for purposes of allocating memory for internal variables, for instance
                %                 % What happens at slave points in terms of stresses is of no concern for the method)
                %                 %                 [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.etaNON,ECMdata.etaNONder,...
                %                 %                     ECMdata.etaNONder2] =...
                %                 %                     DiscreteECM_adaptWEIGHTS(A_internalFORCES_ECM,setPoints,wRED,DATA_interp,OPERFE.wSTs) ;
                %
                %                 DATAoffline= DefaultField(DATAoffline,'Hyperreduction_Separate_Slave_contribution',1);  % Introduced 1-Sept-2025, see
                %                 % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/08_PLAST_adapECM.mlx
                %
                %                 if DATAoffline.Hyperreduction_Separate_Slave_contribution == 1
                %                     % Version before 1-Sept-2025
                %                     [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.DATA_regress_eta_der] =...
                %                         DiscreteECM_adaptWEIGHTSfst(A_internalFORCES_ECM,setPoints,wRED,DATA_interp_ECM,OPERFE.wSTs) ;
                %                 else
                %                     [ECMdata.setPoints, ECMdata.setPoints_slv,ECMdata.wRED,ECMdata.wRED_slv,ECMdata.DATA_regress_eta_der] =...
                %                         DiscreteECM_adaptWEIGHTSfstNS(A_internalFORCES_ECM,setPoints,wRED,DATA_interp_ECM,OPERFE.wSTs) ;
                %                 end
                %
                %
                %                 ECMdata.setElements = large2smallREP(ECMdata.setPoints,DATA.MESH.ngaus) ;
                %                 disp(['Master element(s) =',num2str(ECMdata.setElements(:)')])
                %                 ECMdata.setElements_slv = large2smallREP(ECMdata.setPoints_slv,DATA.MESH.ngaus) ;
                %                 disp(['Slave element(s) =',num2str(ECMdata.setElements_slv(:)')])
            otherwise
                ECMdata.DATA_regress_eta_der =[] ;
        end
        
end

% else
%     error('Option not available (2-July-2025)')
%     % ****************************************************************
%     % Function for computing the position of the integration points
%     % ****************************************************************
%     Nst = OTHER_output.Nst ;
%     NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
%     if  ~exist(NAME_MODES_FOLDER)
%         mkdir(NAME_MODES_FOLDER)
%     end
%     DATA_GENGAUSS.NameFileMesh_ECM = [NAME_MODES_FOLDER,NAME_BASE,'DECMpoints'] ;
%     DATALOC = [] ;
%
%     DATA_GENGAUSS.NameFileMesh_FINT = [NAME_MODES_FOLDER,NAME_BASE,'InternalForceModes'] ;
%     DATALOC = [] ;
%
%     DATA_GENGAUSS.NameFileMesh_CECM = [NAME_MODES_FOLDER,NAME_BASE,'CECMpoints'] ;
%     DATALOC = [] ;
%
%
%
%     [ECMdata] = ContinuousECM(BstRED_l,BasisPone,DATA,OPERFE.wSTs,DATAoffline,DATA_GENGAUSS,...
%         MESH,Nst) ;
% end

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
%disp(['Number of PK2-stress modes = ',num2str(size(BasisStwo,2))])
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
save(NAMEOFFLINEstore,'ECMdata','BasisU','BasisStwo','DATA','DATAoffline','BstRED_reactions' )



diary off

