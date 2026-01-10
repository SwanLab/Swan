function   [MESH,MATPRO_cl,OPERHROM_cl,Fbody_cl,Ftrac_cl,DISP_CONDITIONS_cl,INICOND,DATAHROM_cl,...
    OTHER_data,BasisUall_cl,DistanceCentroid] ...
    =  PreInputDataCLUSTER_HROM(DATAHROM,DATAFE,DIRICHLET,NEUMANN,INITIAL_CONDITIONS,...
    ECMdata_cluster,BasisU_cluster,BasisStwo_cluster,OTHER_INPUTS)
% See
% THIS FUNCTION PREPARES THE DATA NEEDED TO RUN THE FE LARGE STRAIN(HROM)
% ANALYSIS.  It is a modification of its FE counterpart,
% PreProcessInputDataDyn1.m

%
OTHER_data =[] ;
if nargin == 0
    load('tmp.mat')
elseif nargin == 8
    OTHER_INPUTS = []  ;
end


nclusters = length(ECMdata_cluster) ;
MATPRO_cl = cell(1,nclusters) ;
OPERHROM_cl = cell(1,nclusters) ;
Fbody_cl.a = [] ;
Fbody_cl.U = cell(1,nclusters) ; ;

Ftrac_cl.a = [] ;
Ftrac_cl.U = cell(1,nclusters) ;

if isstruct(BasisU_cluster)
    BasisUall_cl.BASIS =  [];
    BasisUall_cl.coeff =  cell(1,nclusters) ;
    
else
    BasisUall_cl = cell(1,nclusters) ;
end

DATAHROM_cl.COMMON = DATAHROM ;
DATAHROM_cl.VAR = cell(1,nclusters) ; ;

DISP_CONDITIONS_cl =  cell(1,nclusters) ;
load(DATAFE.FE_VARIABLES_NAMEstore,'MATPRO','OTHER_output','MESH','OPERFE') ;
MATPRO_FE = MATPRO ;
%DATAHROM_INI = DATAHROM ;

OTHER_INPUTS = DefaultField(OTHER_INPUTS,'TransClustDATA',[]) ;
OTHER_INPUTS.TransClustDATA = DefaultField(OTHER_INPUTS.TransClustDATA,'DistanceCentroid',[]) ;
if isempty(OTHER_INPUTS.TransClustDATA.DistanceCentroid)
    
    DistanceCentroid.C2 = zeros(nclusters,1) ;
    DistanceCentroid.Ct_BasisU = cell(nclusters,nclusters) ;
else
    DistanceCentroid = OTHER_INPUTS.TransClustDATA.DistanceCentroid ;
end

for icluster = 1:nclusters
    %MATPRO = MATPRO_FE ;
    %   DATAHROM = DATAHROM_INI ;
    
%     if icluster == 553
%         disp('Borrar esto ')
%     end
    
    DATAHROM_loc =[] ;
    ECMdata = ECMdata_cluster{icluster} ;
    if ~isstruct(BasisU_cluster)
        BasisU = BasisU_cluster{icluster} ;
    end
    
    
    
    
    
    
    
    % ---------------
    % RECOVER MATPRO
    %---------------
    MATPRO = GetMaterialPropertiesHROM(ECMdata,DATAFE,MATPRO_FE,DATAHROM,OTHER_output)  ;
    
    
    
    
    %%% BOUNDARY CONDITIONS
    % ------------------------
    if icluster == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RB_MOTION',[]) ;
        OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RVE_like_boundary_conditions',0) ;
        if OTHER_INPUTS.RVE_like_boundary_conditions == 1
            % Boundary with conditions defined by the disp. gradient
            % DIRICHLET.MACRODEF
            [DOFrFE,dRfe] = DirichletCONDtime_homogZERO(DIRICHLET,DATAHROM_cl.COMMON,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ; ;
        else
            % Standard
            [DOFrFE,dRfe] = DirichletCONDtime(DIRICHLET,DATAHROM_INI,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ; ;
        end
        
        DOFlFE = 1:DATAFE.MESH.ndof ;
        DOFlFE(DOFrFE) = []  ;
        OTHER_data.DOFlFE = DOFlFE ;
        OTHER_data.DOFrFE = DOFrFE ;
        dRfe = DefaultField(dRfe,'RIGID_BODY_MOTION',[]) ;
        DISP_CONDITIONS = [] ;
        DISP_CONDITIONS.RIGID_BODY_MOTION = dRfe.RIGID_BODY_MOTION ;
        dRfe.RIGID_BODY_MOTION = [] ;
        
        %%%%%%%% DEFINITION OF BasisUall
        % --------------------------------------
        BasisU_r = dRfe.U;   % Constrained part
        
        
        nmodesUr = size(BasisU_r,2) ;
    end
    
    
    
    
    % Precomputed data CENTROIDS
    % --------------------------
    
    
    if isempty(OTHER_INPUTS.TransClustDATA.DistanceCentroid)
        if ~isempty(OTHER_INPUTS.CENTROIDS_cl)
            DistanceCentroid.C2(icluster)  = sum(OTHER_INPUTS.CENTROIDS_cl(DOFlFE,icluster).^2) ;
            for   jcluster = 1:nclusters
                DistanceCentroid.Ct_BasisU{jcluster,icluster}  =   OTHER_INPUTS.CENTROIDS_cl(DOFlFE,jcluster)'*BasisU ;
            end
        end
        
    end
    
    
    if icluster == 1
        if isstruct(BasisU_cluster)
            nmodesDOFl_all = size(BasisU_cluster.BASIS,2) ;
            nmodesDOFr = size(BasisU_r,2) ;
            nmodesALL = nmodesDOFl_all + nmodesDOFr  ;
            BasisUall_cl.BASIS = zeros(DATAFE.MESH.ndof,nmodesALL)  ;
            ROWS = DOFlFE ; 
            COLS = 1:size(BasisU_cluster.BASIS,2)  ;
            BasisUall_cl.BASIS(ROWS,COLS) = BasisU_cluster.BASIS ;
            ROWS = DOFrFE ;
            COLS = COLS(end)+1:nmodesALL ;
            BasisUall_cl.BASIS(ROWS,COLS) = BasisU_r ;
        end
    end
    
    
    if isstruct(BasisU_cluster)
        nmodesUl = size(BasisU_cluster.coeff{icluster},2) ;
        nmodesUall = nmodesUr + nmodesUl ;
        BasisUall_cl.coeff{icluster} = zeros(nmodesALL,nmodesUall) ;
        ROWS = 1:nmodesDOFl_all ;  COLS = 1:nmodesUl ;
        BasisUall_cl.coeff{icluster}(ROWS,COLS) =  BasisU_cluster.coeff{icluster} ;
        ROWS = ROWS(end)+1:nmodesALL ;  COLS = COLS(end)+1:nmodesUall ;
        BasisUall_cl.coeff{icluster}(ROWS,COLS) = eye(nmodesUr,nmodesUr) ;
        BasisUall = [] ;
    else
        nmodesUl = size(BasisU,2) ;
        nmodesUall = nmodesUr + nmodesUl ;
        BasisUall = zeros(DATAFE.MESH.ndof,nmodesUall) ;
        COLS  = 1:nmodesUl ;
        BasisUall(DOFlFE,COLS) = BasisU ;
        COLS  = nmodesUl+1:nmodesUall ;
        BasisUall(DOFrFE,COLS) = BasisU_r ;
    end
    
    DATAHROM_loc.MESH =  [] ;
    DATAHROM_loc.MESH.ndof = nmodesUall;
    % ------------------------------------------
    % DIRICHLET CONDITIONS REDUCED-ORDER MODEL
    % -------------------------------------------------
    DISP_CONDITIONS.DOFl = 1:nmodesUl ;
    DISP_CONDITIONS.DOFr = nmodesUl+1:nmodesUall;
    dR.U = eye(nmodesUr);
    dR.a = dRfe.a ;
    DISP_CONDITIONS.dR = dR  ;
    %%% OPERATORS
    % -----------
    
    nF = DATAFE.MESH.ndim^2 ;
    if isempty(ECMdata.setPoints)
        % Interpolation  (  we are using the Continuous ECM )
        % --------------------------------------------------------------------
        error('Option not available yet ')
        %    OPERHROM.Bst =  InterpolationGaussVariablesECM(BstRED,ECMdata,DATAFE.MESH.ngaus_STRESS,nF) ;
        %    setIndices = small2large(setPointsElement,nF) ;
    else
        setIndices = small2large(ECMdata.setPoints,nF) ;
        %   OPERHROM.Bst = BstRED(setIndices,:) ;
    end
    
    if isstruct(BasisU_cluster)
        if icluster == 1
            BstRED_ini = OPERFE.Bst*BasisUall_cl.BASIS ;
            OPERHROM.Bst = BstRED_ini(setIndices,:)*BasisUall_cl.coeff{icluster} ;
        else
            OPERHROM.Bst = BstRED_ini(setIndices,:)*BasisUall_cl.coeff{icluster} ;
        end
    else
        OPERHROM.Bst = OPERFE.Bst(setIndices,:)*BasisUall ;
    end
    
    OPERHROM.wSTs =  ECMdata.wRED;
    % OPERHROM.IDENTITY_F = [] ; % OPERFE.IDENTITY_F(setIndices);
    DATAHROM_loc.MESH.nstrain = DATAFE.MESH.nstrain ;
    DATAHROM_loc.MESH.ngaus_STRESS = DATAFE.MESH.ngaus_STRESS;
    
    DATAHROM_loc.MESH.ngausT =  length(OPERHROM.wSTs) ; % Total number of Gauss points
    DATAHROM_loc.MESH.ndofSTRESS =   DATAHROM_loc.MESH.ngausT*DATAHROM_loc.MESH.nstrain ; % Total number of Gauss points
    DATAHROM_loc.MESH.ndof = nmodesUall;
    
    % % 5.  MASS MATRIX
    % % -----------------------
    % densGLO = repmat(MATPRO.dens',ngaus_RHS,1) ;
    % densGLO = densGLO(:) ;
    % densGLO_W = CompWeightDiag(densGLO.*wSTs_RHS,ndim) ;
    % NstW = densGLO_W*Nst ;
    % OPERFE.M = NstW'*Nst ;
    
    % --------------------------------------------------
    % 6. Nodal body forces due to gravity (SELF-WEIGHT)
    % --------------------------------------------------
    %Fbody.U =  BasisUall'*OTHER_output.Fbody_U ;
    
    if icluster==1
        FbodyFE.U =  OTHER_output.Fbody_U ;
        FbodyFE.a = ones(size(DATAHROM_cl.COMMON.STEPS)) ;    % Space-Time decomposition.  Gravity is the same for all time steps
        
        if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
            % Rotation of body forces
            %     % -----------------------
            FbodyFE =    RotateForces(Fbody,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATAHROM_cl.COMMON,MESH) ;
        end
        Fbody_cl.a = FbodyFE.a ;
        if isstruct(BasisU_cluster)
            FbodyFE.U = BasisUall_cl.BASIS'*FbodyFE.U ;
        end
    end
    
    if isstruct(BasisU_cluster)
        Fbody_cl.U{icluster}  =  BasisUall_cl.coeff{icluster}'*FbodyFE.U  ;
    else
        Fbody_cl.U{icluster} =   BasisUall'*FbodyFE.U  ;
    end
    
    
    
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 8. Neumann (natural) boundary conditions :   NODAL FORCES due to distributed loads
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISTRIBUTED LOADS  and point loads (NOTE: these ARE NOT follower loads)
    % ------------------------
    if icluster == 1
        FtracFE = NeumannCONDtime(NEUMANN,DATAHROM_cl.COMMON,DATAFE.MESH.ndim,MESH,OTHER_output.NstT_W_N_boundaries) ;
        if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
            % Rotation of body forces
            % -----------------------
            FtracFE =    RotateForces(Ftrac,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATAHROM_cl.COMMON,MESH) ;
        end
        Ftrac_cl.a = FtracFE.a ;
        
        if isstruct(BasisU_cluster)
            FtracFE.U = BasisUall_cl.BASIS'*FtracFE.U ;
        end
    end
    
    if isstruct(BasisU_cluster)
        Ftrac_cl.U{icluster}  =  BasisUall_cl.coeff{icluster}'*FtracFE.U  ;
    else
        Ftrac_cl.U{icluster} =   BasisUall'*FtracFE.U  ;
    end
    
    
    
    
    
    %%%%
    MATPRO_cl{icluster} = MATPRO ;
    OPERHROM_cl{icluster} = OPERHROM ;
    %  Fbody_cl{icluster} = Fbody ;
    %  Ftrac_cl{icluster} = Ftrac ;
    if  ~isempty(BasisUall)
        BasisUall_cl{icluster} = BasisUall ;
    end
    DISP_CONDITIONS_cl{icluster} =DISP_CONDITIONS ;
    % BasisUr_cl{icluster} = BasisU_r ;
    
    
    DATAHROM_cl.VAR{icluster} =DATAHROM_loc ;
end

% *************************
% Compressing BasisU_r
% *************************
% TOL_BLOCK = zeros(length(BasisUr_cl),1)' ;
% DATAsvd=[];
% [U,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(BasisUr_cl,TOL_BLOCK,DATAsvd) ;
% BasisU_r_cl = [] ; BasisU_r_cl.U = U ;
% coeff = bsxfun(@times,V',S) ;
% BasisU_r_cl.coeff = coeff ;
OTHER_data.BasisU_r    = BasisU_r ;

% -------------------------
% 9. INITIAL CONDITIONS
% --------------------------
INICOND = [] ;


% ndof = prod(size(MESH.COOR)) ;
% DOFl = setdiff(1:ndof,DOFr) ;
% INICOND.DISP = zeros(ndof,1) ;
% INICOND.VELOC = zeros(ndof,1) ;
%
% INICOND.DISP(DOFr) = dR.U*dR.a(:,1) ;  % Compatibility with   boundary conditions
% if INITIAL_CONDITIONS.DISP.ISZERO == 0
%     error('You must specify the .mat file with the input initial displacement ')
% end
%
% if INITIAL_CONDITIONS.VELOC.ISZERO == 0
%     error('You must specify the .mat file with the input initial velocity ')
% end


