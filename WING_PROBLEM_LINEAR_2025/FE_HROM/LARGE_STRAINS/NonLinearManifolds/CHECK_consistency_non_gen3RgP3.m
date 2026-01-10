function CHECK_consistency_non_gen3RgP(NameParamStudyHROM,NAME_BASE,DATALOC,DATALOC1,other_INPUTS)
%==========================================================================
% Adaptation of CHECK_consistency_non_gen3Rg to elastostatic problems
% JAHO, 12-Oct-2025, Sunday, Balmes 185, Barcelona
%
% CHECK_consistency_non_gen3Rg
%
% PURPOSE
%   End–to–end **consistency check** for a nonlinear manifold HROM that may
%   also include **reaction-force** outputs and (optionally) **homogenized
%   PK1 stresses** from boundary reactions. The routine compares HROM
%   reconstructions against the corresponding FE (FOM) snapshot fields that
%   were used during training, using the affine-BC manifold setting:
%
%       d  = (A Φ) · τ(qL)  +  D^U · qR .
%
%   It reports relative errors, plots diagnostic curves, and (if enabled)
%   evaluates reaction histories and their induced homogenized stresses.
%
% WHAT IT CHECKS / PRODUCES
%   • **Displacements**: reconstruct HROM displacements via τ(qL) and the
%     affine basis (AΦ, D^U); compare to FE snapshot displacements.
%   • **PK2 stresses**: compare HROM vs FE snapshot PK2 stress fields
%     (projected SVD reconstructions).
%   • **Reactions** (optional): build HROM reaction outputs using the
%     reduced reaction operator BstRED_reactions and (possibly) **q-dependent
%     ECM weights**; compare to FE reactions on a selected boundary face/DOF.
%   • **Homogenized PK1 stress** (optional): compute from boundary reactions
%     (both FE and HROM) using a geometric assembly matrix over boundary nodes.
%   • **Diagnostics**: per-project relative errors and quick plots (norms at
%     ECM points, reaction histories, homogenized stress components).
%   • **Report**: a short textual summary printed to CONSISTENTY.txt.
%
% HIGH-LEVEL WORKFLOW
%   1) Load OFFLINE data (./SNAPSHOTS/NAME_BASE/OFFLINE.mat):
%        BasisU, BasisU_r{proj}, ECMdata, DATAoffline, OTHER_output,
%        BstRED_reactions (if available).
%      Extract DISP_CONDITIONS_HROM and the **precompiled τ evaluator**
%      (OTHER_output.DATAHROM.DATA_evaluateTAU_and_DER).
%   2) For each project case:
%        a) Load FE snapshot SVD reconstructions (U,S,V) for
%           DISPLACEMENTS, PK2, PK1, and RESID (reactions).
%        b) Load HROM snapshot SVD reconstructions for the same fields.
%        c) Reconstruct HROM displacements via:
%               [qL; qR]  →  [τ(qL); qR]  →  d = (AΦ)τ(qL) + D^U qR
%           using ReconsDispl_MHROM and project-specific BasisU_r{iproj}.
%        d) Determine **ECM weights**:
%               – fixed: wRED (vector), or
%               – parametric: wRED(q) via ECMdata.wRED.DATA_regress.
%        e) Compute relative Frobenius errors in displacements
%           (and, if desired, PK2 stresses).
%        f) If requested (DATALOC1.COMPARE_STRESS_ECM_POINTS), plot stress
%           norms at ECM points for FE vs HROM.
%        g) If a boundary face/DOF is specified (DATALOC1.FACE_TO_ANALIZE,
%           DATALOC1.DOF_TO_ANALIZE), compare **reaction histories**:
%               FE:   sum over selected face DOFs from RESID
%               HROM: Output_HROM_depend_reactions(DATA, BstRED_reactions,
%                                                  wECM, PK1stress)
%        h) If requested (DATALOC1.ComputeHomogeneizedStresses==1), assemble
%           **homogenized PK1 stress** from boundary reactions using a geometry
%           matrix X_homog_stress and divide by the domain volume.
%
% KEY INPUTS
%   NameParamStudyHROM : Base name for HROM results in ./SNAPSHOTShrom/<name>_param_*
%   NAME_BASE          : Base name for OFFLINE + FE snapshots in ./SNAPSHOTS/<name>_param_*
%   DATALOC            : Struct with .Trajectories() (returns INPUTS_PARAMETERS)
%   DATALOC1           : Optional controls:
%                         • FACE_TO_ANALIZE (integer face id)
%                         • DOF_TO_ANALIZE  (component on that face, 1..ndim)
%                         • COMPARE_STRESS_ECM_POINTS (true/false)
%                         • ComputeHomogeneizedStresses (true/false)
%   other_INPUTS       : Optional; may override FE snapshots base via
%                        .NAME_BASE_FE_SNAPSHOTS
%
% IMPORTANT FILES / STRUCTS EXPECTED
%   • ./SNAPSHOTS/NAME_BASE/OFFLINE.mat  →  {ECMdata, BasisU, BasisU_r, DATAoffline,
%                                            OTHER_output, BstRED_reactions}
%   • ./SNAPSHOTS/<NAME>_param_k/INFO_SNAPSHOTS.mat   (FE)
%   • ./SNAPSHOTShrom/<HROM>_param_k/INFO_SNAPSHOTS.mat  (HROM)
%   • Each snapshot MAT includes SVD factors for DISP, PK2STRESS, PK1STRESS,
%     and RESID (reactions) used here for fast reconstruction.
%
% FORMULAS / CONVENTIONS
%   • Affine Dirichlet reconstruction (slides):
%         d = (AΦ) τ(qL) + D^U qR .
%   • Relative error for a field X:
%         err = ||X_HROM − X_FE||_F / ||X_FE||_F .
%   • Homogenized PK1 (2D example in the code):
%         P̄ = (1/|Ω|) · X_homog_stress · R_bnd ,
%     where R_bnd stacks boundary reaction DOFs, and X_homog_stress depends on
%     boundary-node coordinates.
%
% DIAGNOSTICS / OUTPUT
%   • Prints parameters and errors to CONSISTENTY.txt.
%   • Stores per-project errors in TableInfoProj.{DISP,PK2}.
%   • Generates simple plots for stress norms at ECM points, reaction histories,
%     and homogenized PK1 components (when requested).
%
% PITFALLS / TIPS
%   • Ensure FE/HROM directories correspond to the **same** projects and step
%     counts; otherwise comparisons will be meaningless.
%   • BasisU_r is **project dependent** (cell array); always use BasisU_r{iproj}.
%   • If ECM weights are parametric (struct), evaluate them at the same qL used
%     for τ(qL) before computing reaction-dependent quantities.
%   • Units: homogenized stress requires dividing by the domain volume; check
%     consistency of length/force units (MPa vs N/mm², etc.).
%
% HISTORY
%   • 2024-09-11  Fast B-spline evaluation path.
%   • 2025-07-28  Input reshuffle & robustness tweaks.
%   • 2025-10-02  Adds **reaction-force** and **homogenized PK1** checks.
%   • 2025-10-08  Comments updated to reflect affine-BC manifold reconstruction.
%==========================================================================


%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
    close all
    
end
% if nargin == 3
%     other_INPUTS = [] ;
% end
%------------------------------------------------------------
delete('CONSISTENTY.txt')
diary 'CONSISTENTY.txt'
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
%
%-------------------------
% FE info                *
% ------------------------
other_INPUTS = DefaultField(other_INPUTS,'NAME_BASE_FE_SNAPSHOTS',NAME_BASE) ;
NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,other_INPUTS.NAME_BASE_FE_SNAPSHOTS];



NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];

% HROM info
NAMEsnap_baseHROM = [cd,filesep,'SNAPSHOTShrom',filesep,NameParamStudyHROM,'_param_'];

DATAloc= feval(DATALOC.Trajectories) ;

CASES = 1:length(DATAloc.INPUTS_PARAMETERS) ;  % Number of training projects


% Let us construct first the matrix the snapshots of each project
load(NAMEOFFLINEstore,'ECMdata','BasisU','DATAoffline','OTHER_output','BasisU_r','BstRED_reactions','Kll')

DISP_CONDITIONS_HROM = OTHER_output.DISP_CONDITIONS;

disp('REDUCED-ORDER DATA')
disp('*****************************')
if isfield(DATAoffline,'errorDISP')
    disp(['  DISP svd tolerance (block) =',num2str(DATAoffline.errorDISP)])
end
%disp(['  STRESS svd tolerance (block) =',num2str(DATAoffline.errorSTRESS)])
disp([' FINT svd tolerance   =',num2str(DATAoffline.errorFINT)])
%disp(['ECM  tolerance  =',num2str(DATAoffline.errorECM)])
disp(['Number of disp. modes (extended) = ',num2str(size(BasisU,2))])
if ~isstruct(ECMdata.wRED)
    disp(['Number of ECM points = ',num2str(length(ECMdata.wRED))])
else
    disp(['Subspace Adaptive WEights ECM'])
    disp(['Number of ECM points =',num2str(length(ECMdata.setPoints))]);
    disp(['ECM elements =',num2str(ECMdata.setElements(:)')])
    disp(['ECM approx. function =',ECMdata.wRED.DATA_regress.nameFunctionEvaluate])
end
disp('Nonlinear mappings')
DATA_evaluateTAU_and_DER = OTHER_output.DATAHROM.DATA_evaluateTAU_and_DER;
disp(DATA_evaluateTAU_and_DER.nameFunctionEvaluate)
disp('Number of generalized coordinates')
disp(['Unconstrained = ',num2str(length(DISP_CONDITIONS_HROM.DOFl))])
disp(['Constrained = ',num2str(length(DISP_CONDITIONS_HROM.DOFr))])


TableInfoProj.PK2 = zeros(size(CASES));
TableInfoProj.DISP =zeros(size(CASES));







for iproj = 1:length(CASES)
    % FE
    % -----------------
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_FOLDER_HROM = [NAMEsnap_baseHROM,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    if iproj == 1
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'MATPRO','MESH','OPERFE','DISP_CONDITIONS') ;
        DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
        SMALL_STRAIN_KINEMATICS =  DATA.SMALL_STRAIN_KINEMATICS ;
        
        
        % INFO for REACTIONS
        
        DATALOC1 = DefaultField(DATALOC1,'FACE_TO_ANALIZE',[]) ;
        
        if ~isempty(DATALOC1.FACE_TO_ANALIZE)
            FACE_TO_ANALIZE = DATALOC1.FACE_TO_ANALIZE ;
            NODES_FACES= MESH.NODES_FACES{FACE_TO_ANALIZE} ;
            ndim = size(MESH.COOR,2) ;
            DOFs_face = small2large(NODES_FACES,ndim) ;
            DOFS_selected = DOFs_face(DATALOC1.DOF_TO_ANALIZE:ndim:end) ;
        else
            FACE_TO_ANALIZE = [] ;
            
        end
        
        
        DATALOC1 = DefaultField(DATALOC1,'ComputeHomogeneizedStresses',0) ;
        
        if DATALOC1.ComputeHomogeneizedStresses == 1
            
            AllNodesBoundary = unique(cell2mat(MESH.NODES_FACES(:))) ; %{FACE_TO_ANALIZE} ;
            Xbnd = MESH.COOR(AllNodesBoundary,:) ;
            % Now we define a Mxnstress (nstress = 4 for  2D) containing
            % rows that depends on Xbnd
            if size(MESH.COOR,2) == 2
                nSTRESS = 4 ;
                X_homog_stress  = zeros(nSTRESS,2*size(Xbnd,1)) ;
                % PK1 homog. stress will be calculated as the sum of  size(Xbdn,1) matrices of the form
                % P_node_i = [React_1*X_1
                % React_2*X_2
                % React_1*X_2
                % React_2*X_1 ]
                %  = [X1   0
                %     0    X_2
                %     X2    0
                %     0     X1]    *
                % [React_1
                %  React_2]
                X_homog_stress(1,1:2:end) =    Xbnd(:,1) ;
                X_homog_stress(2,2:2:end) =    Xbnd(:,2) ;
                X_homog_stress(3,1:2:end) =    Xbnd(:,2) ;
                X_homog_stress(4,2:2:end) =    Xbnd(:,1) ;
                
            else
                error('Option not implemented yet')
            end
        end
    end
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    STRESS_LOC = cell(1,length(NAME_SNAP_loc)) ;
    REACTIONS_FE =cell(1,length(NAME_SNAP_loc)) ;
    PK1STRESS_FE = cell(1,length(NAME_SNAP_loc)) ;
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        DISP_LOC{iloc} = DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
        STRESS_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
        STRESS_LOC{iloc} = STRESS_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        if SMALL_STRAIN_KINEMATICS == 0
            
            PK1STRESS_FE{iloc} = bsxfun(@times,SNAP_cluster.PK1STRESS.U',SNAP_cluster.PK1STRESS.S)' ;
            PK1STRESS_FE{iloc} = PK1STRESS_FE{iloc}*SNAP_cluster.PK1STRESS.V' ;
            
        else
            PK1STRESS_FE{iloc} =  STRESS_LOC{iloc} ;
        end
        
        REACTIONS_FE{iloc} = bsxfun(@times,SNAP_cluster.RESID.U',SNAP_cluster.RESID.S)' ;
        REACTIONS_FE{iloc} =  REACTIONS_FE{iloc}*SNAP_cluster.RESID.V' ;
        
        
    end
    DISP_LOC_FE= cell2mat(DISP_LOC);
    STRESS_LOC_FE= cell2mat(STRESS_LOC);
    REACTIONS_FE= cell2mat(REACTIONS_FE);
    PK1STRESS_FE =cell2mat(PK1STRESS_FE) ;
    % HROM
    % ----
    NAME_INFO_HROM = [NAME_FOLDER_HROM,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO_HROM,'INFO_SNAPSHOTS')
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    STRESS_LOC = cell(1,length(NAME_SNAP_loc)) ;
    PK1STRESS = cell(1,length(NAME_SNAP_loc)) ;
    
    %  REACTIONS_HROM =cell(1,length(NAME_SNAP_loc)) ;
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        DISP_LOC{iloc} = DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
        
        
        
        STRESS_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
        STRESS_LOC{iloc} = STRESS_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        if SMALL_STRAIN_KINEMATICS == 0
            PK1STRESS{iloc} = bsxfun(@times,SNAP_cluster.PK1STRESS.U',SNAP_cluster.PK1STRESS.S)' ;
            PK1STRESS{iloc} = PK1STRESS{iloc}*SNAP_cluster.PK1STRESS.V' ;
        else
            PK1STRESS{iloc} =  STRESS_LOC{iloc} ;
        end
        
        %     REACTIONS_HROM{iloc} = bsxfun(@times,SNAP_cluster.RESID.U',SNAP_cluster.RESID.S)' ;
        %     REACTIONS_HROM{iloc} =  REACTIONS_HROM{iloc}*SNAP_cluster.RESID.V' ;
        
        
    end
    
    
    % REACTIONS_HROM= cell2mat(REACTIONS_HROM);
    PK1STRESS= cell2mat(PK1STRESS);
    
    DISP_LOC_HROM= cell2mat(DISP_LOC);
    
    STRESS_LOC_HROM= cell2mat(STRESS_LOC);
    
    [DISP_LOC_HROM,qL] = ReconsDispl_MHROM(DISP_LOC_HROM,DISP_CONDITIONS_HROM,DATA_evaluateTAU_and_DER,BasisU,DISP_CONDITIONS,BasisU_r{iproj}) ;
    
    
    % Evaluate integration weights
    if isstruct(ECMdata.wRED)
        qW = qL(ECMdata.wRED.IndexDOFl_q,:) ;
        WeigthsECM = feval(ECMdata.wRED.DATA_regress.nameFunctionEvaluate,qW,ECMdata.wRED.DATA_regress) ;
    else
        WeigthsECM = repmat(ECMdata.wRED,1,size(qL,2)) ;
    end
    
    %  DISP_LOC_HROM = BasisUall*DISP_LOC_HROM_proj ;
    % Reconstruction stresses
    %   STRESS_LOC_HROM = RECONS_PK2STRESS.coeff*STRESS_LOC_HROM ;
    %   STRESS_LOC_HROM = RECONS_PK2STRESS.BASIS*STRESS_LOC_HROM ;
    % ERROR displacements
    % --------------------
    nERROR_disp = norm(DISP_LOC_HROM-DISP_LOC_FE,'fro') ;
    nDISP =  norm(DISP_LOC_FE,'fro') ;
    errorLOC = nERROR_disp/nDISP ;
    TableInfoProj.DISP(iproj) =  errorLOC ;
    disp('************************************************')
    disp(['Project = ',num2str(iproj)])
    disp('************************************************')
    disp(['REL. DISP. ERROR  = ',num2str(errorLOC)])
    
    
    %     disp('Examining the error at snapshot level')
    %     DISP_LOC_FE_proj = BasisUall'*DISP_LOC_FE ;
    %     DOFl_HROM_ext = 1:size(qL_extend,1) ;
    %     qL_extend_FE = DISP_LOC_FE_proj(DOFl_HROM_ext,:) ;
    %
    %     n_qL_extend = sqrt(sum(qL_extend.^2,1)) ;
    %     n_qL_extend_FE = sqrt(sum(qL_extend_FE.^2,1)) ;
    %     figure(56+iproj)
    %     hold on
    %     xlabel('Snapshot')
    %     ylabel('Norm qL extendend')
    %     title(['PROJECT = ',num2str(iproj)])
    %     h1 = plot(n_qL_extend,'DisplayName','HROM','LineStyle','--') ;
    %     h2 = plot(n_qL_extend_FE,'DisplayName','FE') ;
    %     legend show;
    
    %  DATALOC1.COMPARE_STRESS_ECM_POINTS = false ;
    DATALOC1 = DefaultField(DATALOC1,'COMPARE_STRESS_ECM_POINTS',true) ;
    
    if DATALOC1.COMPARE_STRESS_ECM_POINTS
        
        IndexexECM_to_plot_stresses = 1:length(ECMdata.setPoints);
        
        
        Index_ECMpoint = ECMdata.setPoints(IndexexECM_to_plot_stresses) ;
        nstrain = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.MESH.nstrain;
        indexes_STRESS_FE = small2large(Index_ECMpoint,nstrain) ;
        indexes_STRESS_HROM = small2large(IndexexECM_to_plot_stresses,nstrain) ;
        
        
        stressPOINT_HROM = sqrt(sum(STRESS_LOC_HROM(indexes_STRESS_HROM,:).^2,1)) ;
        stressPOINT_FE = sqrt(sum(STRESS_LOC_FE(indexes_STRESS_FE,:).^2,1)) ;
        LL = 'Norm Stress (MPa)' ;
        
        
        figure(845+iproj)
        hold on
        xlabel('Snapshot')
        ylabel(LL)
        if  length(Index_ECMpoint) == 1
            title(['Stress point ',num2str(Index_ECMpoint)])
        else
            title(['Norm stresses all ECM points ','PROJECT = ',num2str(iproj)])
        end
        h1 = plot(stressPOINT_HROM,'DisplayName','HROM','LineStyle','--') ;
        h2 = plot(stressPOINT_FE,'DisplayName','FE') ;
        legend show;
        
    end
    
    %
    %     nERROR_stress = norm(STRESS_LOC_HROM-STRESS_LOC_FE,'fro') ;
    %     nSTRESS =  norm(STRESS_LOC_FE,'fro') ;
    %     errorLOC = nERROR_stress/nSTRESS ;
    
    %  disp(['REL. PK2 STRESS. ERROR  = ',num2str(errorLOC)]) ;
    
    %
    TableInfoProj.PK2(iproj) =  errorLOC ;
    
    if ~isempty(FACE_TO_ANALIZE)
        
        figure(211+iproj)
        hold on
        title(['PROJ =',num2str(iproj),' Reaction FACE ',num2str(DATALOC1.FACE_TO_ANALIZE),' direction ',num2str(DATALOC1.DOF_TO_ANALIZE)])
        xlabel('Step')
        ylabel('REaction (MN)')
        REACTIONS_time = sum(REACTIONS_FE(DOFS_selected,:),1) ;
        plot(REACTIONS_time,'DisplayName','FE')
        
        % COMPUTING REACTIONS, HROM
        REACTIONS_time_HROM = Output_HROM_depend_reactions(DATA,BstRED_reactions,WeigthsECM,PK1STRESS) ;
        
        
        plot(REACTIONS_time_HROM,'DisplayName','HROM')
        legend show
        
    end
    
    
    if DATALOC1.ComputeHomogeneizedStresses == 1 && ~isempty(BstRED_reactions)
        
        figure(665+iproj)
        hold on
        title(['PROJ =',num2str(iproj),' Homogeneized PK1 stress ' ])
        xlabel('Step')
        ylabel('Stress (MPa)')
        DOFsBND = small2large(AllNodesBoundary,DATA.MESH.ndim) ;
        
        PK1stress_homg =  X_homog_stress*REACTIONS_FE(DOFsBND,:)  ;
        VOL = sum(OPERFE.wSTs) ;
        PK1stress_homg = PK1stress_homg/VOL ;
        
        % HROM
        % -----------------------------
        PK1stress_homg_HROM = Output_HROM_depend_reactions(DATA,BstRED_reactions,WeigthsECM,PK1STRESS) ;
        
        
        for  iSTRESS = 1:size(PK1stress_homg,1)
            
            
            plot(PK1stress_homg(iSTRESS,:),'DisplayName',['P_',num2str(iSTRESS),  ' FE'])
            plot(PK1stress_homg_HROM(iSTRESS,:),'DisplayName',['P_',num2str(iSTRESS),  ' HROM'])
        end
        
        legend show
        
    end
    
    
    figure(347+iproj)
    hold on
    title(['Project = ',num2str(iproj),' Latent variables '])
    xlabel('Snapshot')
    ylabel('q')
    
    for  iLATENT = 1:size(qL,1)
    plot(qL(iLATENT,:),'DisplayName',['HROM q_',num2str(iLATENT)],'LineStyle','--') ;
  %  h2 = plot(qL(2,:),'DisplayName','Plastic q','LineStyle','--') ;
    
    
    qL_FE = BasisU(:,iLATENT)'*Kll*DISP_LOC_FE(DISP_CONDITIONS.DOFl,:) ;
     plot(qL_FE,'DisplayName',['FE q_',num2str(iLATENT)],'LineStyle','--') ;
    
    end
    
    % FE

    legend show;
    
    
    
    
end

%
% TableInfoProj.PK2 = ERROR_STRESSES ;
% TableInfoProj.DISP = ERROR_DISP ;
% TableInfoProj.TRANSITION = MAX_TRANS_ERROR ;


%TableModes(TableInfoProj) ;

diary off