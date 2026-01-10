function   [MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATA,OTHER_data] ...
    =  PreProcessInputDataDyn1(DATA,PROPMAT,DIRICHLET,NEUMANN,INITIAL_CONDITIONS,OTHER_INPUTS)
%%
% =========================================================================
% PreProcessInputDataDyn1 — Build full FE pack (mesh, operators, BCs, loads)
% =========================================================================
% PURPOSE
%   Assemble all inputs required by the large-strain FE solver from high-level
%   descriptors: mesh & geometric operators, material at Gauss points, mass/
%   damping matrices, gravity/body and surface loads, essential BCs, initial
%   conditions, storage/printing policies, and auxiliary metadata.
%
% SIGNATURE
%   [MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATA,OTHER_data] = ...
%       PreProcessInputDataDyn1(DATA, PROPMAT, DIRICHLET, NEUMANN, INITIAL_CONDITIONS, OTHER_INPUTS)
%
% INPUTS
%   DATA                (struct)  Global run config (problem type, flags, gravity, printing, etc.).
%   PROPMAT             (struct)  Material parameters (possibly multiple materials).
%   DIRICHLET           (struct)  Prescribed-displacement boundary data (time-dependent).
%   NEUMANN             (struct)  Traction boundary data (time-dependent, non-follower).
%   INITIAL_CONDITIONS  (struct)  Optional dINI/vINI vectors (defaults to zero).
%   OTHER_INPUTS        (struct)  Optional extras (alphaD/betaD damping, RB_MOTION, etc.).
%
% OUTPUTS
%   MESH            : COOR, CN, boundary connectivity, face props, centroid, volume, inertia, center of mass.
%   MATPRO          : Material fields at Gauss points (elasticity/plastic data, density).
%   OPERFE          : FE operators: Bst (deformation gradient), wSTs, mass M, (optional) geometric mass Mgeo,
%                     optional hydrostatic follower-load operators, potential-energy operator, etc.
%   Fbody, Ftrac    : Space–time separated body and traction loads (gravity & Neumann).
%   DISP_CONDITIONS : DOFr, DOFl, dR (and rigid-body motion info if present).
%   INICOND         : Initial displacement and velocity vectors.
%   DATA            : Updated config (mesh sizes, nstrain, clustering, printing).
%   OTHER_data      : Aux info (geom properties, shape-function storage, damping flags, etc.).
%
% WHAT THE FUNCTION DOES
%   1) Mesh & sizes
%      • Build mesh (GeometryMesh or GeometryMeshTILED); set ndim, nnode, nelem.
%      • Decide strain vector length nstrain (3/4 in 2D, 6 in 3D); allow 4-component plane-strain option.
%   2) Geometric operators
%      • Compute Bst_F, Nst, quadrature weights/points, face properties, centroid, volume/inertia.
%      • Optionally store shape-function matrices for diagnostics.
%   3) Material at Gauss points
%      • Map PROPMAT to Gauss fields via MaterialPropertiesIntegPoints (sets MATPRO, DATA, INICOND).
%   4) Mass (and geometric mass) matrices
%      • Build consistent M; optionally build Mgeo when ComputeGeometricMassMatrix is true.
%      • Compute center of mass and initial potential energy operator for gravity.
%   5) Damping
%      • Rayleigh-like αM + βK: assemble αM here; defer βK to linear-stiffness stage if requested.
%   6) Storage planning
%      • Estimate memory; compute clusters and filenames for snapshot storage and GiD printing.
%   7) Essential BCs (Dirichlet)
%      • Build DOFr/dR via DirichletCONDtime_GENERAL (supports periodicity/RVE BCs and RB motion).
%   8) Loads
%      • Body forces from gravity; rotate by RB motion if present.
%      • Surface tractions from NEUMANN via NeumannCONDtime; rotate if RB motion is present.
%      • Optional follower loads (e.g., hydrostatic) and waterline mesh.
%   9) Initial conditions
%      • Set dINI/vINI (zeros by default); optional consistency checks for dynamics.
%
% KEY DEPENDENCIES
%   GeometryMesh / GeometryMeshTILED, GeometricMatricesFun, MaterialPropertiesIntegPoints,
%   DirichletCONDtime_GENERAL, NeumannCONDtime, RotateForces, HYDROST_computeOPER,
%   ClusterComputeSize, DefaultField, CompWeightDiag.
%
% PRACTICAL NOTES
%   • DATA.NO_USE_Deformation_gradient_in_Small_Strains and DATA.SMALL_STRAIN_KINEMATICS
%     decide the strain measure size used by Bst_F (F-based vs small-strain).
%   • Set DATA.LIMIT_mbytes_matrices to manage snapshot clustering.
%   • For dynamic runs with Newmark-Bossak, βK damping part is flagged for computation after K is available.
%   • Gravity is assumed along a single axis; potential energy is computed if enabled.
%
% USAGE EXAMPLE
%   [MESH,MATPRO,OPERFE,Fbody,Ftrac,DBC,INICOND,DATA,OTHER] = ...
%       PreProcessInputDataDyn1(DATA, PROPMAT, DIRICH, NEUMANN, INICOND, struct('alphaD',0,'betaD',0));
%
% .mlx references:
%   • README_plasticity.mlx — plasticity workflow notes
%   • TowerSimulation.mlx — example on initial displacements due to forces at t=0
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================


%
OTHER_data =[] ;
if nargin == 0
    load('tmp.mat')
end


% ----------------------------
%%%%%%%%%%%%%%%%%
% ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB, Connectivity faces,
% normals ....
DATA = DefaultField(DATA,'TILED_MESH',[]) ; 
if  isempty(DATA.TILED_MESH)
    MESH = GeometryMesh(DATA) ;
else
    MESH = GeometryMeshTILED(DATA) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA = DefaultField(DATA,'StrainStressWith4Components',0) ;
[nnode,ndim ]= size(MESH.COOR) ;% Number of nodes
[nelem,nnodeE ]= size(MESH.CN) ; % Number of elements

ndim = size(MESH.COOR,2)  ;
if ndim==2 && DATA.StrainStressWith4Components == 0
    nstrain = 3;
elseif ndim==2 && DATA.StrainStressWith4Components == 1 && strcmp(DATA.typePROBLEM,'pstrain')
    nstrain =4 ;
elseif ndim==2 && DATA.StrainStressWith4Components == 1
    error('Plane stress is not compatible with the option DATA.StrainStressWith4Components = ')
else
    nstrain = 6 ;
    typePROBLEM ='3D' ;
end
DATA.MESH.nstrain  =nstrain ;

DATA = DefaultField(DATA,'posgp_given',[]) ;
DATA = DefaultField(DATA,'weights_given',[]) ;
MESH.posgp_given = DATA.posgp_given ;
MESH.weights_given = DATA.weights_given ;






% 3. GEOMETRIC MATRICES (AND RENUMBERING)
% Geometric matrices
%-------------------------------------
% DATA = DefaultField(DATA,'RenumberElementsForEficiency',1) ;
% if  DATA.RenumberElementsForEficiency == 1
%     disp('Renumering elements')
%     [~,IndicesRenumberingElements]  = sort(MESH.CN(:,1)) ;
%     MESH.CN = MESH.CN(IndicesRenumberingElements,:) ;
%     MATPRO.celasglo = MATPRO.celasglo(:,:,IndicesRenumberingElements) ;
%     MATPRO.dens = MATPRO.dens(IndicesRenumberingElements) ;
%     MESH.MaterialType = MESH.MaterialType(IndicesRenumberingElements) ;
% else
%     IndicesRenumberingElements = [] ;
% endw
% MESH.IndicesRenumberingElements =IndicesRenumberingElements ;

% 4. FE/HROM OPERATORS,MATRICES
% ---------------------------
DATA = DefaultField(DATA,'NO_USE_Deformation_gradient_in_Small_Strains',1) ;
DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0) ;

if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 0 || DATA.SMALL_STRAIN_KINEMATICS ==0
    
    DATA.MESH.nstrain_F = ndim^2 ;
else
    DATA.MESH.nstrain_F = nstrain ;
    
end



MESH.DATA = DATA;
[Bst_F,wSTs,Nst,wSTs_RHS,NstT_W_N_boundaries,ngaus_RHS,GEOproperties,ngaus_STRESS,IDENTITY_F,posgp,shapef_RHS] ...
    = GeometricMatricesFun(MESH,nstrain) ;
MESH.DATA = [] ;

MESH.PROPERTIES_FACES = GEOproperties.FACES ;
MESH.CENTROID =  GEOproperties.CENTROID ;
MESH.VOLUME =  GEOproperties.VOLUME ;
MESH.INERTIA =  GEOproperties.INERTIA ;

disp(['Total volume = ',num2str(MESH.VOLUME)])

DATA = DefaultField(DATA,'STORE_MATRIX_SHAPE_FUNCTIONS_Nst',0) ;
if DATA.STORE_MATRIX_SHAPE_FUNCTIONS_Nst == 1
    OTHER_data.Nst = Nst ;
end

OTHER_data.GEOproperties = GEOproperties;
DATA.MESH.posgp  =posgp ;
DATA.MESH.ngaus  =ngaus_STRESS ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo, density   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA = DefaultField(DATA,'TYPE_CONSTITUTIVE_MODEL_ALL','SMALL_STRAINS_ELASTIC') ;
DATA.ListFieldInternalVariables = [] ; % Use this variable to specify the list of internal variables
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
INICOND = [] ;

[MATPRO,MESH,DATA,INICOND] = MaterialPropertiesIntegPoints(MESH,DATA,PROPMAT)  ; % 23-March-2023, to be tested


% switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
%     case 'SMALL_STRAINS_ELASTIC'
%         MESH.nstrain=  DATA.MESH.nstrain  ;
%         [MATPRO] = SmallStrainElasticityPROP(MESH,DATA.typePROBLEM,PROPMAT) ;
%         disp('Global elasticity matrix    ...')
%         Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
%         %  Cglo = ConvertBlockDiag(Cglo) ;
%         MATPRO.celasglo = [] ;
%         MATPRO.celasglo = Cglo ;
%
%     case 'NeoHookean'
%         [MATPRO] = NeoHookPROP(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
%         disp('Global elasticity matrix    ...')
%         %     Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
%         %  Cglo = ConvertBlockDiag(Cglo) ;
%         %    MATPRO.celasglo = [] ;
%         %   MATPRO.celasglo = Cglo ;
%     case 'SMALL_STRAINS_J2_PLASTICITY'
%         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
%         [MATPRO,DATA] = SmallStrainJ2PlasticityPROP(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
%         INICOND.YieldStress =  MATPRO.sigmay_0 ;   % Initial condition
%         DATA.ListFieldInternalVariables = {'YieldStress','PlasticStrains','InternalVarStrain'} ;
%         ngausT = size(INICOND.YieldStress,1) ;
%         INICOND.PlasticStrains = zeros(ngausT*DATA.MESH.nstrain,1) ;
%         INICOND.InternalVarStrain = zeros(size(INICOND.YieldStress)) ;
% end



OPERFE.Bst = Bst_F ;   % Matrix such that F = Bst_F*d + ident
OPERFE.IDENTITY_F = IDENTITY_F ;
%OPERFE.wSTs = wSTs ;     % Integration weights at the Gauss points used for evaluating stresses and internal
% forces
% Matrix of integration weights
% _F = ndim^2 ;
% wST = sparse(repmat(wSTs',_F,1)) ;
% wST = wST(:) ;
% wST = spdiags(wST,1,length(wST),length(wST))  ; %
%S = spdiags(Bin,d,m,n) creates an m-by-n sparse matrix S by taking the columns of Bin and placing them along the diagonals specified by d.
OPERFE.wSTs = wSTs ;

DATA.MESH.ngausT =  length(wSTs) ; % Total number of Gauss points
DATA.MESH.ndofSTRESS =   DATA.MESH.ngausT*nstrain; % Total number of Gauss points
DATA.MESH.ngaus_STRESS = ngaus_STRESS;
DATA.MESH.ndof = nnode*ndim;
DATA.MESH.ndim = ndim;

% 5.  MASS MATRIX
% -----------------------
DATA.MESH.ngaus_RHS = ngaus_RHS ;
densGLO = repmat(MATPRO.dens',ngaus_RHS,1) ;
densGLO = densGLO(:) ;
densGLO_W = CompWeightDiag(densGLO.*wSTs_RHS,ndim) ;
NstW = densGLO_W*Nst ;
OPERFE.M = NstW'*Nst ;

DATA = DefaultField(DATA,'ComputeGeometricMassMatrix',false) ; 

if DATA.ComputeGeometricMassMatrix
    densGLO_W = CompWeightDiag(ones(size(densGLO)).*wSTs_RHS,ndim) ;
    NstW = densGLO_W*Nst ;
    OPERFE.Mgeo = NstW'*Nst ;
else
     OPERFE.Mgeo = [] ; 
end




% Mass center (gravity)
CENTER_MASS = zeros(1,ndim) ;
mass_point =  wSTs_RHS.*densGLO;
mass = sum(mass_point) ;
%for idim = 1:3
COOR = MESH.COOR' ;
NstCOOR = Nst*COOR(:) ;
for idim = 1:ndim
    CENTER_MASS(idim) = (mass_point)'*NstCOOR(idim:ndim:end)/mass ;
end

MESH.CENTER_MASS = CENTER_MASS;

% INITIAL POTENTIAL ENERGY (y = 0 )
DATA = DefaultField(DATA,'PotentialEnergyCompute',1) ;
if any(DATA.vGRAVITY)  && DATA.PotentialEnergyCompute == 1
    
    ggg = find(DATA.vGRAVITY ~=0) ;
    if length(ggg) >1
        error('Gravity should act along one coordinate axis')
    end
    idim = ggg;
    P = sum((NstCOOR(idim:ndim:end).*abs(DATA.vGRAVITY(idim))).*mass_point)  ;
    DATA.InitialPotentialEnergy = P ;
    DATA.idim_gravity = idim ;
    
    OPERFE.Nst_potential_energy = abs(DATA.vGRAVITY(idim))*NstW(idim:ndim:end,:) ;
else
    
    OPERFE.Nst_potential_energy = 0 ;
end




%



% 5.A) Damping matrix
% -------------------
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'alphaD',0) ;
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'betaD',0) ;
OTHER_data.ComputeLinearStiffnessMatrixForDamping_beta = 0 ;

 
if OTHER_INPUTS.alphaD == 0
    OPERFE.Ddamp = [] ; %
else
    OPERFE.Ddamp = OTHER_INPUTS.alphaD*OPERFE.M ;
    if OTHER_INPUTS.betaD ~= 0
        disp('The remaining component of the damping matrix will be computed when determing the linear stiffness matrix') ; 
        disp('In  ComputeKiniCheck.m')
        OTHER_data.ComputeLinearStiffnessMatrixForDamping_beta = OTHER_INPUTS.betaD ; 
    end
end



%OPERFE.Ddamp  = alphaBAR*OPERFE.M ;


% COMPUTING NUMBER OF CLUSTERS (FOR STORING PURPOSES)
% ******************************
DATA = DefaultField(DATA,'LIMIT_mbytes_matrices',50) ;
DATA = ClusterComputeSize(DATA.LIMIT_mbytes_matrices,DATA) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RB_MOTION',[]) ;
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'MESHcoarse_FOR_DIRICHLET_BOUND_COND',[]) ;
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RVE_like_boundary_conditions',0) ;
%DATAINPUT = DefaultField(DATAINPUT,'TYPE_BOUNDARY_CONDITIONS','PERIODICITY_ALL_FACES_WITH_CORNERS') ;


% General function (see particular cases inside )
[DOFr,dR,OTHEROUTPUT_DIRICH] = DirichletCONDtime_GENERAL(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
OTHEROUTPUT_DIRICH = DefaultField(OTHEROUTPUT_DIRICH,'DISP_CONDITIONS',[]) ;
%OTHER_data.INFO_PERIODIC_CONDITIONS = OTHEROUTPUT_DIRICH.INFO_PERIODIC_CONDITIONS ;
MESH.INFO_PERIODIC_CONDITIONS =  OTHEROUTPUT_DIRICH.INFO_PERIODIC_CONDITIONS ;

if  isempty(OTHEROUTPUT_DIRICH.DISP_CONDITIONS)
    DISP_CONDITIONS.DOFr = DOFr;
    DOFl = 1:DATA.MESH.ndof ;
    DOFl(DOFr) = []  ;
    DISP_CONDITIONS.DOFl = DOFl ;
    DISP_CONDITIONS.dR = dR;
    dR = DefaultField(dR,'RIGID_BODY_MOTION',[]) ;
    DISP_CONDITIONS.RIGID_BODY_MOTION = dR.RIGID_BODY_MOTION ;
    dR.RIGID_BODY_MOTION = [] ;
else
    DISP_CONDITIONS = OTHEROUTPUT_DIRICH.DISP_CONDITIONS ;
    DISP_CONDITIONS.RIGID_BODY_MOTION  = [] ;
end
% -----------------------------------------------------------------------------------

% 7. Nodal body forces due to gravity (SELF-WEIGHT)

gV = repmat(DATA.vGRAVITY(:),nnode,1) ;
Fbody.U = OPERFE.M*gV ;
Fbody.a = ones(size(DATA.STEPS)) ;
%Fbody.a(1:3) = [0,0.25,0.5];
%disp('borrar lo de arriba')
OTHER_data.Fbody_U = Fbody.U ;

% Fbody.a = ones(size(DATA.STEPS)) ;    % Space-Time decomposition.  Gravity is the same for all time steps

% for idim = 1:ndim
%     gV(idim:ndim:end) = gV(idim:ndim:end).*densGLO;
% end
% Fbody.U = NstW'*gV ;
% %dR = DefaultField(dR,'RIGID_BODY_MOTION',[] ) ;
if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
    % Rotation of body forces
    % -----------------------
    Fbody =    RotateForces(Fbody,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATA,MESH) ;
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Neumann (natural) boundary conditions :   NODAL FORCES due to distributed loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS  and point loads (NOTE: these ARE NOT follower loads)


%
OTHEROUTPUT_DIRICH= DefaultField(OTHEROUTPUT_DIRICH,'Ftrac',[]); %.Ftrac
if isempty(OTHEROUTPUT_DIRICH.Ftrac)
    Ftrac = NeumannCONDtime(NEUMANN,DATA,ndim,MESH,NstT_W_N_boundaries,GEOproperties) ;
else
    % See [DOFr,dR,OTHEROUTPUT_DIRICH] = DirichletCONDtime_GENERAL(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
    Ftrac = OTHEROUTPUT_DIRICH.Ftrac ;
end

if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
    % Rotation of body forces
    % -----------------------
    Ftrac =    RotateForces(Ftrac,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATA,MESH) ;
end

OTHER_data.NstT_W_N_boundaries = NstT_W_N_boundaries ;




% 8. FOLLOWER LOADS (HYDROSTATIC)
% --------------------------------
DATA = DefaultField(DATA,'FOLLOWER_LOADS',[]) ;
OPERFE.HYDRO = [] ;

if ~isempty(DATA.FOLLOWER_LOADS)
    switch  DATA.FOLLOWER_LOADS.TYPE
        case 'HYDROSTATIC'
            
            DATA.FOLLOWER_LOADS.HYDROSTATIC = DefaultField(DATA.FOLLOWER_LOADS.HYDROSTATIC,'INCLUDE_DYNAMIC_PRESSURE',0) ;
            
            % Hydrostatic forces
            [OPERFE.HYDRO,DATA,MESH ]= HYDROST_computeOPER(MESH,DATA)   ;
            %
            DATA.FOLLOWER_LOADS.HYDROSTATIC  = DefaultField(DATA.FOLLOWER_LOADS.HYDROSTATIC,'NAME_MESH_WATERLINE',[]) ;
            if ~isempty(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE)
                % WATERLINE MESH
                % ------------------
                [DATA.WATERLINE_MESH]= ReadMeshFileStr(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE,'READ_MATERIAL_COLUMN',1)  ;
                
            end
            
            
        otherwise
            error('Option not implemented')
    end
    
end



% -------------------------
% 9. INITIAL CONDITIONS
% --------------------------

ndof = prod(size(MESH.COOR)) ;
DOFl = setdiff(1:ndof,DOFr) ;

INITIAL_CONDITIONS = DefaultField(INITIAL_CONDITIONS,'dINI',[]) ;
INITIAL_CONDITIONS = DefaultField(INITIAL_CONDITIONS,'vINI',[]) ;

if isempty(INITIAL_CONDITIONS.dINI)
    INICOND.DISP = zeros(ndof,1) ;
else
    INICOND.DISP = INITIAL_CONDITIONS.dINI ;
end

if isempty(INITIAL_CONDITIONS.vINI)
    INICOND.VELOC = zeros(ndof,1) ;
else
    INICOND.VELOC = INITIAL_CONDITIONS.vINI ;
end


% INITIAL DISPLACEMENT DUE TO GRAVITY
% ONLY WHEN INITIAL_CONDITIONS.dINI IS NOT ZERO, and there are  disp.
% cots
%if all(INITIAL_CONDITIONS.dINI==0) && any(DATA.vGRAVITY)  && DATA.ISDYNAMIC == 1 && ~isempty(DISP_CONDITIONS.DOFr)
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'DATAINPUT',[]) ;
OTHER_INPUTS.DATAINPUT = DefaultField(OTHER_INPUTS.DATAINPUT,'INITIAL_DISPLACEMENTS_CAUSED_BY_FORCES_AT_ZERO',0) ;

InitDispF = OTHER_INPUTS.DATAINPUT.INITIAL_DISPLACEMENTS_CAUSED_BY_FORCES_AT_ZERO ;  % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/03_CYLINDRICAL_TOWER/TowerSimulation.mlx

if all(INITIAL_CONDITIONS.dINI==0)   && DATA.ISDYNAMIC == 1 && ~isempty(DISP_CONDITIONS.DOFr) ...
        && InitDispF ==1   % 25-May-2022 . Include also cases with no gravity, just forces at t = 0
    
    disp(['Computing initial displacement due to applied forces at t = 0'])
    
    INICOND = InitialDisplacementsDueToGravity(DATA,INICOND,DISP_CONDITIONS,Fbody,Ftrac,OPERFE,MATPRO) ;
    OTHER_data.INICOND = INICOND ;
else
    OTHER_data.INICOND = INICOND ;
end






if DATA.ISDYNAMIC == 1
    
    switch DATA.INTEGRATION_SCHEME.TYPE
        case 'NEWMARKbossak'
            a = DATA.INTEGRATION_SCHEME.NEWMARKbossak.aBOSSAK ;
            DATA.TIME_INTEGRATION_PARAMETERS.a = a ;
            if a > 0
                error('aBOSSAK should be less than zero ')
            end
            %  \begin{equation}
            %   \gamma = \dfrac{1}{2} - \aBOSSAK
            %  \end{equation}
            DATA.TIME_INTEGRATION_PARAMETERS.TYPE = 1;
            DATA.TIME_INTEGRATION_PARAMETERS.GAMMA = 0.5-a ;
            %  \begin{equation}
            %   \beta = \dfrac{1}{4} \Par{1 - \aBOSSAK}^2
            %  \end{equation}
            DATA.TIME_INTEGRATION_PARAMETERS.BETA = 0.25*(1-a)^2 ;
        otherwise
            error('Option not implemented')
    end
    
    if ~isempty(dR.U)
        BoundaryConditionsR = dR.U*dR.a(:,1) ;
        InitialConditionsR =  INICOND.DISP(DOFr) ;
        
        AA = norm(BoundaryConditionsR-InitialConditionsR)  ;
        
        if  AA > 1e-14
            disp('Boundary conditions not consistent with initial condition')
            pause
        end
        
    end
    
end
