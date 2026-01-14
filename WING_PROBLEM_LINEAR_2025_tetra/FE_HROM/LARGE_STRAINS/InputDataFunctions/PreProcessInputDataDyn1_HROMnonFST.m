function   [MESH,MATPRO,OPERHROM,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATAHROM,OTHER_data,BasisUall] ...
    =  PreProcessInputDataDyn1_HROMnonFST(DATAHROM,DATAFE,DIRICHLET,NEUMANN,INITIAL_CONDITIONS,ECMdata,BasisU,BasisStwo,...
    OTHER_INPUTS   )
%--------------------------------------------------------------------------
% PreProcessInputDataDyn1_HROMnonFST
% ----------------------------------
% PURPOSE
%   Assemble all inputs required by the ONLINE HROM solver on a τ(q)-manifold:
%   - Enforce Dirichlet conditions (including optional rigid-body motion) and
%     build consistent DOF partitions (DOFl/DOFr) inherited from the OFFLINE run.
%   - Construct the global reduced basis BasisUall compatible with affine BCs.
%   - Build reduced operators (M, optional Ddamp, optional KinternalFORCES_given),
%     ECM-weighted stress operators, and optional hydro/follower-load operators.
%   - Project external actions (body forces, tractions) into the reduced space.
%   - Set initial conditions in reduced coordinates and prepare integration params.
%
% SCOPE / CONTEXT
%   This routine “bridges” OFFLINE artifacts (mesh/operators/materials/bases/ECM)
%   with ONLINE execution. It removes anonymous-evaluator dependencies by
%   installing a precompiled evaluator for manifold ECM master→slave mapping:
%       OPERHROM.DATA_regress_eta_der
%   enabling fast, vectorized reconstruction of slave-point work/stress and
%   their derivatives at runtime.
%
% INPUTS (high level)
%   DATAHROM          : Master configuration for ONLINE run (time grid, printing,
%                       integration scheme/parameters, storage settings, flags).
%   DATAFE            : Full-order model data used during OFFLINE stage (mesh,
%                       FE operators, partitions, GEOproperties, etc.).
%   DIRICHLET         : Callback/data for essential BCs (may include rigid motion).
%   NEUMANN           : Callback/data for natural BCs (distributed/point loads).
%   INITIAL_CONDITIONS: Optional initial d and v (full space) to be projected.
%   ECMdata           : Hyperreduction info (wRED, optional wRED_slv, and
%                       DATA_regress_eta_der for manifold ECM).
%   BasisU, BasisStwo : Reduced bases for displacement and PK2-stress recovery.
%   OTHER_INPUTS      : Options (Rayleigh damping α/β or estimation from ω_n,
%                       RVE-like BCs, hydro settings, etc.).
%
% OUTPUTS
%   MESH, MATPRO          : FE mesh & material data for the run.
%   OPERHROM              : Reduced operators & auxiliaries:
%                           • M           : reduced mass.
%                           • Ddamp       : optional Rayleigh damping.
%                           • Bst (ECM)   : strain/stress operator at selected points.
%                           • wSTs, wRED_slv : ECM weights (masters/slaves).
%                           • DATA_regress_eta_der : precompiled manifold-ECM evaluator.
%                           • KinternalFORCES_given : (linear-only shortcut, optional).
%                           • Nst_potential_energy  : energy operator (optional).
%   Fbody, Ftrac          : Reduced body/traction forces (with optional rotation for RBM).
%   DISP_CONDITIONS       : DOF partitions, affine data, and rigid-motion info.
%   INICOND               : Initial reduced displacement/velocity.
%   DATAHROM (updated)    : Mesh sizes, gauss counts, scheme constants (e.g., Bossak γ, β).
%   OTHER_data            : Aux partitions and bookkeeping (DOFlFE/DOFrFE, reactions data).
%   BasisUall             : Global basis consistent with current affine constraints.
%
% KEY STEPS / WORKFLOW
%   1) Recover FE operators and materials
%      - load(DATAFE.FE_VARIABLES_NAMEstore, 'MATPRO','OTHER_output','MESH','OPERFE')
%      - GetMaterialPropertiesHROM(...) to adjust material/initial data.
%
%   2) Dirichlet boundary conditions & DOF partitions
%      - DirichletCONDtime_GENERAL(...) produces DOFrFE and dRfe (and optionally
%        DISP_CONDITIONS). If DISP_CONDITIONS not provided, enforce consistency
%        with OFFLINE DOFr and build DOFlFE.
%      - Build BasisUall and ECM-weighted operators via:
%          • Def_BasisU_Bst_HROMper (if DISP_CONDITIONS.G is present)
%          • Def_BasisU_Bst_HROMgen (general affine BCs)
%
%   3) Manifold ECM (master/slave mapping)
%      - If ECMdata.DATA_regress_eta_der exists, install it in OPERHROM and store
%        wRED_slv (slave weights). This enables nonlinear slave-point reconstruction
%        from master ECM points without runtime anonymous closures.
%
%   4) Reduced operators and shortcuts
%      - Mass matrix:            OPERHROM.M = Uᵀ M U
%      - Linear shortcut (opt.): OPERHROM.KinternalFORCES_given = Uᵀ K_lin U
%      - Damping:
%          • If DATAHROM.EstimateCoefficientsDampingFromNatFrequencies == 1:
%              EstimateDdampFromNatFrequencies(...)
%          • Else, build Rayleigh damping from OTHER_INPUTS.alphaD/betaD when given.
%
%   5) External actions
%      - Body forces: Fbody.U and time law Fbody.a; rotate if RBM is active; project UᵀF.
%      - Tractions:   Ftrac from NeumannCONDtime(...); rotate if RBM is active; project UᵀF.
%      - Optional potential-energy operator mapped to reduced space if provided by OPERFE.
%
%   6) Hydrostatic follower loads (optional)
%      - If OPERFE.HYDRO exists: HYDRO_HROMoperators(...) to fill OPERHROM.HYDRO.
%
%   7) Initial conditions & integration parameters
%      - Project INITIAL_CONDITIONS.dINI and vINI (if given) with Uᵀ(·); else zeros.
%      - If dynamic and NEWMARKbossak: set a (aBOSSAK), then Γ=0.5−a, Β=¼(1−a)².
%        Optionally, check BC/IC consistency on prescribed DOFs.
%
% FEATURES
%   - Affine-BC-consistent basis construction (BasisUall) preserving training partitions.
%   - ECM restricted operators and optional manifold-ECM with precompiled regressors.
%   - Optional linear internal-force shortcut for purely linear cases.
%   - Optional Rayleigh damping or automatic α/β estimation from target frequencies.
%   - RBM-aware rotation of external actions.
%
% ASSUMPTIONS / CONTRACTS
%   - DOFr used here must match the OFFLINE training DOFr (hard check included).
%   - ECMdata provides consistent sets/weights with those used OFFLINE.
%   - BasisU spans the intended free-DOF subspace after enforcing affine BCs.
%
% COMMON PITFALLS / TIPS
%   - Mismatch in DOFr partitions between OFFLINE and ONLINE will break affine maps.
%   - If enabling manifold-ECM, ensure DATA_regress_eta_der includes all derivatives
%     required by the chosen solver/tangent updates.
%   - When estimating damping from frequencies, verify units and discretization step.
%
% DEPENDENCIES (invoked here)
%   - GetMaterialPropertiesHROM, DirichletCONDtime_GENERAL,
%     Def_BasisU_Bst_HROMgen / Def_BasisU_Bst_HROMper,
%     NeumannCONDtime, RotateForces, EstimateDdampFromNatFrequencies,
%     HYDRO_HROMoperators.
%
% VERSION HISTORY (keep existing references)
%   - Prior comments update: “Comments updated by ChatGPT, 8-Oct-2025”.
%
% AUTHOR
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------%
OTHER_data =[] ;
if nargin == 0
    load('tmp1.mat')
elseif nargin == 8
    OTHER_INPUTS = []  ;
end


% OPERHROM.tauNON = tauNON;   % Nonlinear mapping between amplitudes modes and reduced coordinates
% OPERHROM.tauNONder = tauNONder;  % Derivative of tauNON with respect reduced coordinates
% OPERHROM.tauNONder2 = tauNONder2;  % Derivative of tauNON with respect reduced coordinates


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


% if OTHER_INPUTS.RVE_like_boundary_conditions == 1
%     % Boundary with conditions defined by the disp. gradient
%     % DIRICHLET.MACRODEF
%     [DOFrFE,dRfe] = DirichletCONDtime_homogZERO(DIRICHLET,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ; ;
% else
%     % Standard
%     [DOFrFE,dRfe] = DirichletCONDtime(DIRICHLET,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ; ;
% end

[DOFrFE,dRfe,OTHEROUTPUT_DIRICH] = DirichletCONDtime_GENERAL(DIRICHLET,DATAHROM,DATAFE.MESH.ndim,MESH,OTHER_output.GEOproperties,OTHER_INPUTS) ;
OTHEROUTPUT_DIRICH = DefaultField(OTHEROUTPUT_DIRICH,'DISP_CONDITIONS',[]) ;

if  isempty(OTHEROUTPUT_DIRICH.DISP_CONDITIONS) 
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
    
else
    DISP_CONDITIONS = OTHEROUTPUT_DIRICH.DISP_CONDITIONS ;
    DISP_CONDITIONS.RIGID_BODY_MOTION  = [] ;
    
    OTHER_data.DOFlFE = DISP_CONDITIONS.DOFl ;
    OTHER_data.DOFrFE =  DISP_CONDITIONS.DOFr ;
    dRfe.RIGID_BODY_MOTION = [] ;
end


if isfield(DISP_CONDITIONS,'G')
    [DISP_CONDITIONS,OPERHROM,OTHER_data,BasisUall] =...
        Def_BasisU_Bst_HROMper(BasisU,DATAFE,OTHER_INPUTS,DATAHROM,OTHER_data,OPERFE,ECMdata,DISP_CONDITIONS) ;
else
    % Affine boundary conditions in the FE model
    [DISP_CONDITIONS,OPERHROM,OTHER_data,BasisUall] =...
        Def_BasisU_Bst_HROMgen(DOFrFE,dRfe,BasisU,DATAFE,OTHER_INPUTS,DATAHROM,OTHER_data,OPERFE,ECMdata,DOFlFE,DISP_CONDITIONS) ;
end





OPERHROM.wSTs =  ECMdata.wRED;
ECMdata = DefaultField(ECMdata,'DATA_regress_eta_der',[]) ; 
if ~isempty(ECMdata.DATA_regress_eta_der)
    % This is nonlinear version for the ECM points
    % The idea is that the internal work density of some "slave" points
    % is determined as a function (ECMdata.etaNON) of the internal work density of the ECM
    % points 
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/03_ECMmanifold.mlx
    OPERHROM.DATA_regress_eta_der = ECMdata.DATA_regress_eta_der ; 
%     OPERHROM.etaNONder = ECMdata.etaNONder ; 
%     OPERHROM.etaNONder2 = ECMdata.etaNONder2 ; 
    OPERHROM.wRED_slv   = ECMdata.wRED_slv ; 
else
    OPERHROM.DATA_regress_eta_der = [] ; 
end

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
DATAHROM.MESH.ndof =  length(DISP_CONDITIONS.DOFl) +  length(DISP_CONDITIONS.DOFr); 

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

DISP_CONDITIONS = DefaultField(DISP_CONDITIONS,'RIGID_BODY_MOTION',[]) ; 
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
