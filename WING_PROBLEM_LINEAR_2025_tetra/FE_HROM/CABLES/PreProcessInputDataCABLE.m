function   [MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATA,OTHER_data] ...
    =  PreProcessInputDataCABLE(DATA,PROPMAT,DIRICHLET,NEUMANN,INITIAL_CONDITIONS,OTHER_INPUTS)
% This is an adaptation of
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/PreProcessInputDataDyn1.m
% for cable problems
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/...
%  04_MOORING_PROBLEMS/01_STATICgrav.mlx
% JAHO, 22-June-2022

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
MESH = DATA.MESH1D ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nnode,ndim ]= size(MESH.COOR) ;% Number of nodes
[nelem,nnodeE ]= size(MESH.CN) ; % Number of elements
ndim = size(MESH.COOR,2)  ;
nstrain = 1;
DATA.MESH.nstrain  =nstrain ;

DATA = DefaultField(DATA,'posgp_given',[]) ;
DATA = DefaultField(DATA,'weights_given',[]) ;
MESH.posgp_given = DATA.posgp_given ;
MESH.weights_given = DATA.weights_given ;



% 3. GEOMETRIC MATRICES (AND RENUMBERING)
% Geometric matrices
%-------------------------------------

MESH.DATA = DATA;
[Bst_F,wSTs,Nst,wSTs_RHS,NstT_W_N_boundaries,ngaus_RHS,GEOproperties,ngaus_STRESS,IDENTITY_F,posgp,shapef_RHS] ...
    = GeometricMatricesCABLE(MESH,nstrain) ;
MESH.DATA = [] ;
%
% MESH.PROPERTIES_FACES = GEOproperties.FACES ;
% MESH.CENTROID =  GEOproperties.CENTROID ;
% MESH.VOLUME =  GEOproperties.VOLUME ;
% MESH.INERTIA =  GEOproperties.INERTIA ;
%
% disp(['Total volume = ',num2str(MESH.VOLUME)])

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
%switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
%     case 'SMALL_STRAINS_ELASTIC'
%         MESH.nstrain=  DATA.MESH.nstrain  ;
%         [MATPRO] = SmallStrainElasticityPROPcable(MESH,PROPMAT) ;
%         disp('Global elasticity matrix    ...')
%      %   Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
%         Cglo = zeros(size(MATPRO.celasglo,1)*DATA.MESH.ngaus,1) ;
%         for igaus = 1:DATA.MESH.ngaus
%             Cglo(igaus:DATA.MESH.ngaus:end) = MATPRO.celasglo ;
%         end
%         %  Cglo = ConvertBlockDiag(Cglo) ;
%         MATPRO.celasglo = [] ;
%         MATPRO.celasglo = Cglo ;
%    case 'NONLINEAR_SMALL'
[MATPRO] = SmallStrainNonLinearCABLE(MESH,PROPMAT) ;
%    MATPRO = PROPMAT ;
%end



OPERFE.Bst = Bst_F ;   % Matrix such that F = Bst_F*d + ident
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

% CROSS-SECTIONAL AREA
AREA = repmat(MATPRO.AREA',ngaus_STRESS,1) ;
AREA = AREA(:) ;
MATPRO.AREA = AREA ;

% % Mass center (gravity)
CENTER_MASS = zeros(1,ndim) ;
mass_point =  wSTs_RHS.*densGLO;
mass = sum(mass_point) ;
% %for idim = 1:3
COOR = MESH.COOR' ;
NstCOOR = Nst*COOR(:) ;
for idim = 1:ndim
    CENTER_MASS(idim) = (mass_point)'*NstCOOR(idim:ndim:end)/mass ;
end
%
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

if OTHER_INPUTS.alphaD == 0
    OPERFE.Ddamp = [] ; %
else
    OPERFE.Ddamp = OTHER_INPUTS.alphaD*OPERFE.M ;
    if OTHER_INPUTS.betaD ~= 0
        error('Option not implemented')
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
% OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RB_MOTION',[]) ;
%
% %DIRICHLET = DefaultField(DIRICHLET,'MACRODEF',[]) ;
% OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RVE_like_boundary_conditions',0) ;
% if OTHER_INPUTS.RVE_like_boundary_conditions == 1
%     % Boundary with conditions defined by the disp. gradient
%     % DIRICHLET.MACRODEF
%     [DOFr,dR] = DirichletCONDtime_homogZERO(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
% else
% Standard
[DOFr,dR] = DirichletCONDtimeCABLE(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
%end

DATA = DefaultField(DATA,'DO_NOT_APPLY_DIRICHLET_CONDITION_FIRST_STEP',0) ;
if DATA.DO_NOT_APPLY_DIRICHLET_CONDITION_FIRST_STEP == 1
    for itraj = 1:size(dR.a,1)
        istepZER0=  2;
        dR.a(:,istepZER0) = 0 ;
    end
end




DISP_CONDITIONS.DOFr = DOFr;
DOFl = 1:DATA.MESH.ndof ;
DOFl(DOFr) = []  ;
DISP_CONDITIONS.DOFl = DOFl ;
DISP_CONDITIONS.dR = dR;
% dR = DefaultField(dR,'RIGID_BODY_MOTION',[]) ;
% DISP_CONDITIONS.RIGID_BODY_MOTION = dR.RIGID_BODY_MOTION ;
% dR.RIGID_BODY_MOTION = [] ;


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
% if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
%     % Rotation of body forces
%     % -----------------------
%     Fbody =    RotateForces(Fbody,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATA,MESH) ;
% end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Neumann (natural) boundary conditions :   NODAL FORCES due to distributed loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS  and point loads (NOTE: these ARE NOT follower loads)
% ------------------------
Ftrac = Fbody ;
for iloads = 1:length(Ftrac)
    Ftrac(iloads).a = zeros(size( Ftrac(iloads).a)) ;
end



%; %NeumannCONDtime(NEUMANN,DATA,ndim,MESH,NstT_W_N_boundaries,GEOproperties) ;

% if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
%     % Rotation of body forces
%     % -----------------------
%     Ftrac =    RotateForces(Ftrac,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATA,MESH) ;
% end
%
% OTHER_data.NstT_W_N_boundaries = NstT_W_N_boundaries ;



%
% % 8. FOLLOWER LOADS (HYDROSTATIC)
% % --------------------------------
% DATA = DefaultField(DATA,'FOLLOWER_LOADS',[]) ;
%  OPERFE.HYDRO = [] ;
%
%  if ~isempty(DATA.FOLLOWER_LOADS)
%      switch  DATA.FOLLOWER_LOADS.TYPE
%          case 'HYDROSTATIC'
%
%              DATA.FOLLOWER_LOADS.HYDROSTATIC = DefaultField(DATA.FOLLOWER_LOADS.HYDROSTATIC,'INCLUDE_DYNAMIC_PRESSURE',0) ;
%
%              % Hydrostatic forces
%              [OPERFE.HYDRO,DATA,MESH ]= HYDROST_computeOPER(MESH,DATA)   ;
%              %
%              DATA.FOLLOWER_LOADS.HYDROSTATIC  = DefaultField(DATA.FOLLOWER_LOADS.HYDROSTATIC,'NAME_MESH_WATERLINE',[]) ;
%              if ~isempty(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE)
%                  % WATERLINE MESH
%                  % ------------------
%                  [DATA.WATERLINE_MESH]= ReadMeshFileStr(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE,'READ_MATERIAL_COLUMN',1)  ;
%
%              end
%
%
%          otherwise
%              error('Option not implemented')
%      end
%
%  end



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
% OTHER_INPUTS.DATAINPUT = DefaultField(OTHER_INPUTS.DATAINPUT,'INITIAL_DISPLACEMENTS_CAUSED_BY_FORCES_AT_ZERO',DATA.ISDYNAMIC) ;
%
% InitDispF = OTHER_INPUTS.DATAINPUT.INITIAL_DISPLACEMENTS_CAUSED_BY_FORCES_AT_ZERO ;  % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/03_CYLINDRICAL_TOWER/TowerSimulation.mlx
%
% if all(INITIAL_CONDITIONS.dINI==0)   && DATA.ISDYNAMIC == 1 && ~isempty(DISP_CONDITIONS.DOFr) ...
%         && InitDispF ==1   % 25-May-2022 . Include also cases with no gravity, just forces at t = 0
%
%     disp(['Computing initial displacement due to applied forces at t = 0'])
%
%     INICOND = InitialDisplacementsDueToGravity(DATA,INICOND,DISP_CONDITIONS,Fbody,Ftrac,OPERFE,MATPRO) ;
%     OTHER_data.INICOND = INICOND ;
% else
OTHER_data.INICOND = INICOND ;
% end






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
        
        if  AA > 1e-10
            disp(['Difference between prescribed conditions and initial conditions equal to  ',num2str(AA),'; check that the number of elements is the same in both cases.  PRESS ENTER to go on (or abort)'])
            pause
        end
        
    end
    
end


%
%
% INICOND.VELOC = zeros(ndof,1) ;
%
%
% if INITIAL_CONDITIONS.DISP.ISZERO == 0
%     error('You must specify the .mat file with the input initial displacement ')
% end
%
% if INITIAL_CONDITIONS.VELOC.ISZERO == 0
%     error('You must specify the .mat file with the input initial velocity ')
% end
%
%
