function   [MESH,MATPRO,OPERHROM,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATAHROM,OTHER_data,BasisUall] ...
    =  PreProcessInputDataDyn1_HROM(DATAHROM,DATAFE,DIRICHLET,NEUMANN,INITIAL_CONDITIONS,ECMdata,BasisU,BasisStwo,OTHER_INPUTS)
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


% ---------------
% RECOVER MATPRO
%---------------
load(DATAFE.FE_VARIABLES_NAMEstore,'MATPRO','OTHER_output','MESH','OPERFE') ;
% Getting material properties for the HROM
[MATPRO,INICOND] = GetMaterialPropertiesHROM(ECMdata,DATAFE,MATPRO,DATAHROM,OTHER_output)  ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RB_MOTION',[]) ;

%[DOFrFE,dRfe] = DirichletCONDtime(DIRICHLET,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ;

OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RVE_like_boundary_conditions',0) ;
if OTHER_INPUTS.RVE_like_boundary_conditions == 1
    % Boundary with conditions defined by the disp. gradient
    % DIRICHLET.MACRODEF
    [DOFrFE,dRfe] = DirichletCONDtime_homogZERO(DIRICHLET,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ; ;
else
    % Standard
    [DOFrFE,dRfe] = DirichletCONDtime(DIRICHLET,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ; ;
end




III =  setdiff(DOFrFE,DATAFE.DOFr) ;
if ~isempty(III)
    error('The constrained DOFs should be the same that the one used for training')
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
BasisU_r = full(dRfe.U);   % Constrained part
nmodesUr = size(BasisU_r,2) ;
nmodesUl = size(BasisU,2) ;
nmodesUall = nmodesUr + nmodesUl ;
BasisUall = zeros(DATAFE.MESH.ndof,nmodesUall) ;
COLS  = 1:nmodesUl ;

OTHER_INPUTS = DefaultField(OTHER_INPUTS,'BasisU_allDOFS',[]); 

if ~isempty(OTHER_INPUTS.BasisU_allDOFS)
    % Situations in which the DOFr for training and testing are not the
    % same
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/02_HROMstand.mlx
    BasisU = SVDT(OTHER_INPUTS.BasisU_allDOFS(DOFlFE,:)) ; 
end



BasisUall(DOFlFE,COLS) = BasisU ;
%DATA_BasisUAll.DOFlFE = DOFlFE ;
%DATA_BasisUAll.COLS_DOFlFE = DOFlFE ;
%------------
COLS  = nmodesUl+1:nmodesUall ;
BasisUall(DOFrFE,COLS) = BasisU_r ;
DATAHROM.MESH.ndof = nmodesUall;
%------------
%DATA_BasisUAll.DOFrFE = DOFrFE ;
%DATA_BasisUAll.COLS_DOFrFE = COLS ;
%DATA_BasisUAll.BasisU_r = BasisU_r;
OTHER_data.BasisU_r = BasisU_r;

% ------------------------------------------
% DIRICHLET CONDITIONS REDUCED-ORDER MODEL
% -------------------------------------------------OTHER_data
DISP_CONDITIONS.DOFl = 1:nmodesUl ;
DISP_CONDITIONS.DOFr = nmodesUl+1:nmodesUall;
dR.U = eye(nmodesUr);
dR.a = dRfe.a ;
DISP_CONDITIONS.dR = dR  ;
%%% OPERATORS
% -----------
BstRED = OPERFE.Bst*BasisUall ;
DATAFE = DefaultField(DATAFE,'SMALL_STRAIN_KINEMATICS',0) ; %  =

DATAFE = DefaultField(DATAFE,'NO_USE_Deformation_gradient_in_Small_Strains',1) ; %  =

if DATAFE.SMALL_STRAIN_KINEMATICS == 1 && DATAFE.NO_USE_Deformation_gradient_in_Small_Strains == 1
    nF = DATAFE.MESH.nstrain ;
else
    nF = DATAFE.MESH.ndim^2 ;
end


if isempty(ECMdata.setPoints)
    % Interpolation  (  we are using the Continuous ECM )
    % --------------------------------------------------------------------
    OPERHROM.Bst =  InterpolationGaussVariablesECM(BstRED,ECMdata,DATAFE.MESH.ngaus_STRESS,nF) ;
    setIndices = small2large(1:length(ECMdata.setElements),nF) ;  % This operation should be written more efficiently !!! 3-DEc-2021
    OPERHROM.IDENTITY_F =  OPERFE.IDENTITY_F(setIndices);
else
    setIndices = small2large(ECMdata.setPoints,nF) ;
    OPERHROM.Bst = BstRED(setIndices,:) ;
    if  ~isempty(OPERFE.IDENTITY_F)
        OPERHROM.IDENTITY_F =  OPERFE.IDENTITY_F(setIndices); % This operation should be removed !!! 3-DEc-2021
    else
        OPERHROM.IDENTITY_F = [] ;
    end
end
OPERHROM.wSTs =  ECMdata.wRED;

if  size(OPERHROM.wSTs,2) > 1
    error('THIS METHOD HAS PROVED UNRELIABLE...IT IS UNSTABLE')
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/02_ELASTODYNAMICS.mlx
    % THIS IS THE TAILORED HROM  (DIFFERENT WEIGHT FOR EACH MODE)
    % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/01_VIABILITY.mlx
    % jaho, 3-jAN-2024
    % ---------------------------------------
    OPERHROM.BstW =  OPERHROM.Bst ;
    for imode = 1:size(OPERHROM.wSTs,2)
        WlocPerMode = OPERHROM.wSTs(:,imode) ;
        for istrain = 1:nF
            OPERHROM.BstW(istrain:nF:end,imode) = OPERHROM.Bst(istrain:nF:end,imode).*WlocPerMode ;
        end
    end
    
end

DATAHROM.MESH.nstrain = DATAFE.MESH.nstrain ;
DATAHROM.MESH.ngaus_STRESS = DATAFE.MESH.ngaus_STRESS;

DATAHROM.MESH.ngausT =  length(OPERHROM.wSTs) ; % Total number of Gauss points
DATAHROM.MESH.ndofSTRESS =   DATAHROM.MESH.ngausT*DATAHROM.MESH.nstrain ; % Total number of Gauss points
DATAHROM.MESH.ndof = size(BasisUall,2);

% % 5.  MASS MATRIX
% % -----------------------
OPERHROM.M = BasisUall'*OPERFE.M*BasisUall ;
OTHER_output = DefaultField(OTHER_output,'DATA',[]) ;
OTHER_output.DATA = DefaultField(OTHER_output.DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff',0) ;


if OTHER_output.DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==1
    % This is for purely linear problems
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/03_CYLINDRICAL_TOWER/AcceleratingCodeLinear.mlx
    OPERHROM.KinternalFORCES_given  =  BasisUall'*OPERFE.KinternalFORCES_given*BasisUall ;
else
    OPERHROM.KinternalFORCES_given = [] ;
end



% DAMPING MATRIX
% 5.A) Damping matrix
% -------------------

DATAHROM = DefaultField(DATAHROM,'EstimateCoefficientsDampingFromNatFrequencies',0) ; % 17-May-2024
%OTHER_INPUTS = DefaultField(OTHER_INPUTS,'DampingRatioAllFrequencies',0.25) ; % 17-May-2024
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/02_HROMstand.mlx


if DATAHROM.EstimateCoefficientsDampingFromNatFrequencies ==1
    [OPERHROM.Ddamp,OTHER_data.alphaD,OTHER_data.betaD]...
        = EstimateDdampFromNatFrequencies(OPERHROM,OTHER_output,BasisUall,DISP_CONDITIONS,DATAHROM) ; 
else
    OTHER_INPUTS = DefaultField(OTHER_INPUTS,'alphaD',0) ;
    OTHER_INPUTS = DefaultField(OTHER_INPUTS,'betaD',0) ;
    
    if OTHER_INPUTS.alphaD == 0
        OPERHROM.Ddamp = [] ; %
    else
        OPERHROM.Ddamp = OTHER_INPUTS.alphaD*OPERHROM.M ;
        if OTHER_INPUTS.betaD ~= 0
            error('Option not implemented')
        end
    end
    
end





% --------------------------------------------------
% 6. Nodal body forces due to gravity (SELF-WEIGHT)
% --------------------------------------------------
%Fbody.U =  BasisUall'*OTHER_output.Fbody_U ;

Fbody.U =  OTHER_output.Fbody_U ;
Fbody.a = ones(size(DATAHROM.STEPS)) ;    % Space-Time decomposition.  Gravity is the same for all time steps

OTHER_data.FbodyDOFr.a = Fbody.a  ;       % For computing reactions
OTHER_data.FbodyDOFr.U = Fbody.U  ;

%Fbody.a(1:3) = [0,0.25,0.5];
%disp('borrar lo de arriba')


if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
    % Rotation of body forces
    %     % -----------------------
    Fbody =    RotateForces(Fbody,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATAHROM,MESH) ;
end
Fbody.U =   BasisUall'*Fbody.U  ;



%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Neumann (natural) boundary conditions :   NODAL FORCES due to distributed loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS  and point loads (NOTE: these ARE NOT follower loads)
% ------------------------
Ftrac = NeumannCONDtime(NEUMANN,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.NstT_W_N_boundaries) ;

OTHER_data.FtracDOFr.a = Ftrac.a  ;       % For computing reactions
OTHER_data.FtracDOFr.U = Ftrac.U  ;



if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
    % Rotation of body forces
    % -----------------------
    Ftrac =    RotateForces(Ftrac,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATAHROM,MESH) ;
end


Ftrac.U = BasisUall'*Ftrac.U ;


% FOR COMPUTING POTENTIAL ENERGY  (with respect to x(idim) = 0)

OPERFE = DefaultField(OPERFE,'Nst_potential_energy',[]) ;
if ~isempty(OPERFE.Nst_potential_energy)
    OPERHROM.Nst_potential_energy = OPERFE.Nst_potential_energy*BasisUall ;
end



%   FOLLOWER LOADS (HYDROSTATIC)
% --------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/HYDROSTATICS/HYDROST_computeOPER.m
OPERFE = DefaultField(OPERFE,'HYDRO',[]) ;
if ~isempty(OPERFE.HYDRO)
    [OPERHROM.HYDRO,DATAHROM] = HYDRO_HROMoperators(OPERFE.HYDRO,DATAFE,DATAHROM,OTHER_INPUTS.ECMdata_press,BasisUall) ;
end



% -------------------------
% 9. INITIAL CONDITIONS
% --------------------------
% OTHER_INPUTS = DefaultField(OTHER_INPUTS,'DATAINPUT',[]) ;
% OTHER_INPUTS.DATAINPUT = DefaultField(OTHER_INPUTS.DATAINPUT,'INITIAL_DISPLACEMENTS_CAUSED_BY_FORCES_AT_ZERO',DATA.ISDYNAMIC) ;
%
% if
%

INITIAL_CONDITIONS = DefaultField(INITIAL_CONDITIONS,'dINI',[]) ;
INITIAL_CONDITIONS = DefaultField(INITIAL_CONDITIONS,'vINI',[]) ;








if isempty(INITIAL_CONDITIONS.dINI)
    INICOND.DISP = zeros(DATAHROM.MESH.ndof,1) ;
else
    INICOND.DISP = BasisUall'*INITIAL_CONDITIONS.dINI ;
end

if isempty(INITIAL_CONDITIONS.vINI)
    INICOND.VELOC = zeros(DATAHROM.MESH.ndof,1) ;
else
    INICOND.VELOC = BasisUall'*INITIAL_CONDITIONS.vINI ;
end
if DATAHROM.ISDYNAMIC == 1
    
    switch DATAHROM.INTEGRATION_SCHEME.TYPE
        case 'NEWMARKbossak'
            a = DATAHROM.INTEGRATION_SCHEME.NEWMARKbossak.aBOSSAK ;
            DATAHROM.TIME_INTEGRATION_PARAMETERS.a = a ;
            if a > 0
                error('aBOSSAK should be less than zero ')
            end
            %  \begin{equation}
            %   \gamma = \dfrac{1}{2} - \aBOSSAK
            %  \end{equation}
            DATAHROM.TIME_INTEGRATION_PARAMETERS.TYPE = 1;
            DATAHROM.TIME_INTEGRATION_PARAMETERS.GAMMA = 0.5-a ;
            %  \begin{equation}
            %   \beta = \dfrac{1}{4} \Par{1 - \aBOSSAK}^2
            %  \end{equation}
            DATAHROM.TIME_INTEGRATION_PARAMETERS.BETA = 0.25*(1-a)^2 ;
        otherwise
            error('Option not implemented')
    end
    
    if  ~isempty(dR.U)
        BoundaryConditionsR = dR.U*dR.a(:,1) ;
        InitialConditionsR =  INICOND.DISP(DISP_CONDITIONS.DOFr) ;
        
        AA = norm(BoundaryConditionsR-InitialConditionsR)  ;
        
        if  AA > 1e-10
            disp('Boundary conditions not consistent with initial condition')
            pause
        end
    end
    
end
