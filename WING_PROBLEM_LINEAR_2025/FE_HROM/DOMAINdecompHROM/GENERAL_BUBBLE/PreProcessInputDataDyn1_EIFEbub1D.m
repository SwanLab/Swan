function   [MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATA,OTHER_data] ...
    =  PreProcessInputDataDyn1_EIFEbub1D(DATA,PROPMAT,DIRICHLET,NEUMANN,INITIAL_CONDITIONS,OTHER_INPUTS)
% Adaptation of PreProcessInputDataDyn1_EIFEbub for 1D problems
% JAHO, 13-May-2024, Honest Greens, Tuset, Barcelona
% --------------------------------------
OTHER_data =[] ;
if nargin == 0
    load('tmp2.mat')
end


% ----------------------------
%%%%%%%%%%%%%%%%%
% ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB, Connectivity faces,
% normals ....
DATA.RenumberElementsForEficiency = 0 ;  % THIS WAY THE NUMBERING OF THE PRE-PROCESS FILES IS EQUAL TO THE POST-PROCESS FILES
% HOWEVER, THE CODE WILL PROVE LESS EFFICIENT (NOT BANDED B MATRICES).
% JAHO, 1-APRIL-2023
MESH =  DATA.MESH_GIVEN ; % GeometryMeshUNC(DATA,OTHER_INPUTS) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA = DefaultField(DATA,'StrainStressWith4Components',0) ;
[nnode,ndim ]= size(MESH.COOR) ;% Number of nodes
[nelem,nnodeE ]= size(MESH.CN) ; % Number of elements

nstrain = PROPMAT(1).EIFE_prop.MESH.nstrain ; % nstrain
ndimFINE= size(PROPMAT(1).EIFE_prop.MESH.COOR,2) ; % nstrain

DATA.MESH.nstrain  =nstrain ;

DATA = DefaultField(DATA,'NO_USE_Deformation_gradient_in_Small_Strains',1) ;

if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 0 || DATA.SMALL_STRAIN_KINEMATICS ==0
    
    DATA.MESH.nstrain_F = ndimFINE^2 ;
else
    DATA.MESH.nstrain_F = nstrain ;
    
end
DATA.MESH.ndimFINE = ndimFINE ;

%-----------------------------------------------------------------------------

[Bst_F,wSTs,Nst,wSTs_RHS,ngaus_RHS,ngaus_STRESS,IDENTITY_F,posgp,MESH,DATA,KstiffLINEAR] ...
    = GeometricMatricesFun_EIFEbub1D(MESH,nstrain,PROPMAT,DATA) ;



MESH.DATA = [] ;

% MESH.PROPERTIES_FACES = GEOproperties.FACES ;
% MESH.CENTROID =  GEOproperties.CENTROID ;
% MESH.VOLUME =  GEOproperties.VOLUME ;
% MESH.INERTIA =  GEOproperties.INERTIA ;

% disp(['Total volume = ',num2str(MESH.VOLUME)])

DATA = DefaultField(DATA,'STORE_MATRIX_SHAPE_FUNCTIONS_Nst',0) ;
if DATA.STORE_MATRIX_SHAPE_FUNCTIONS_Nst == 1
    OTHER_data.Nst = Nst ;
end

%OTHER_data.GEOproperties = GEOproperties;
DATA.MESH.posgp  =posgp ;
DATA.MESH.ngaus_STRESS  =ngaus_STRESS ;  % Assumed the same for all elements
DATA.MESH.ngaus_RHS = ngaus_RHS ;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo, density   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA = DefaultField(DATA,'TYPE_CONSTITUTIVE_MODEL_ALL','SMALL_STRAINS_ELASTIC') ;
DATA.ListFieldInternalVariables = [] ; % Use this variable to specify the list of internal variables
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
INICOND = [] ;
switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
    case 'SMALL_STRAINS_ELASTIC'
        MESH.nstrain=  DATA.MESH.nstrain  ;
        [MATPRO] = SmallStrainElasticityPROP_EIFE(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
        
        
    case 'NeoHookean'
        error('Option not implemented')
        %   [MATPRO] = NeoHookPROP(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
        %  disp('Global elasticity matrix    ...')
        %     Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
        %  Cglo = ConvertBlockDiag(Cglo) ;
        %    MATPRO.celasglo = [] ;
        %   MATPRO.celasglo = Cglo ;
    case 'SMALL_STRAINS_J2_PLASTICITY'
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
        %  &     error('Option not implemented yet')
        [MATPRO,DATA] = SmallStrainJ2PlasticityPROP_EIFE(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
        INICOND.YieldStress =  MATPRO.sigmay_0 ;   % Initial condition
        DATA.ListFieldInternalVariables = {'YieldStress','PlasticStrains','InternalVarStrain'} ;
        ngausT = size(INICOND.YieldStress,1) ;
        INICOND.PlasticStrains = zeros(ngausT*DATA.MESH.nstrain,1) ;
        INICOND.InternalVarStrain = zeros(size(INICOND.YieldStress)) ;
end



OPERFE.Bst = Bst_F ;   % Matrix such that F = Bst_F*d + ident
OPERFE.KstiffLINEAR = KstiffLINEAR ;   %  See
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

OPERFE.IDENTITY_F = IDENTITY_F ;

OPERFE.wSTs = wSTs ;

DATA.MESH.ngausT =  length(wSTs) ; % Total number of Gauss points
DATA.MESH.ndofSTRESS =   DATA.MESH.ngausT*nstrain; % Total number of Gauss points
DATA.MESH.ngaus_STRESS = ngaus_STRESS;
DATA.MESH.ndof = length(DATA.MESHextended.DOFS_TO_KEEP);
DATA.MESH.ndim = []; % This is a ambiguos variable now
DATA.MESH.ndimSP = ndim; % This is a ambiguos variable now
DATA.MESH.NDOFS_pernode_max = (DATA.MESHextended.NDOFS_pernode); % This is a ambiguos variable now


% 5.  MASS MATRIX
% -----------------------
%DATA.MESH.ngaus_RHS = ngaus_RHS ;
ndimFINE = size(posgp,1) ;
densGLO = MATPRO.dens ;
densGLO_W = CompWeightDiag(densGLO.*wSTs_RHS,ndimFINE) ;
NstW = densGLO_W*Nst ;
OPERFE.M = NstW'*Nst ;

% % ----------------------------------------------------
% % Mass center (gravity)
% %----------------------------------------------------
% CENTER_MASS = zeros(1,ndim) ;
% mass_point =  wSTs_RHS.*densGLO;
% mass = sum(mass_point) ;
% %for idim = 1:3
% COOR = MESH.COOR' ;
% % NstCOOR = Nst(:,DATA.MESHextended.DOFS_bLOC)*COOR(:) ;
% % for idim = 1:ndim
% %     CENTER_MASS(idim) = (mass_point)'*NstCOOR(idim:ndim:end)/mass ;
% % end
%
% MESH.CENTER_MASS = CENTER_MASS;

% INITIAL POTENTIAL ENERGY (y = 0 )

OPERFE.Nst_potential_energy = 0 ;




%



% 5.A) Damping matrix
% -------------------
% OTHER_INPUTS = DefaultField(OTHER_INPUTS,'alphaD',0) ;
% OTHER_INPUTS = DefaultField(OTHER_INPUTS,'betaD',0) ;
%
% if OTHER_INPUTS.alphaD == 0
%     OPERFE.Ddamp = [] ; %
% else
%     OPERFE.Ddamp = OTHER_INPUTS.alphaD*OPERFE.M ;
%     if OTHER_INPUTS.betaD ~= 0
%         error('Option not implemented')
%     end
 










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

% General function (see particular cases inside )
disp(['Dirichlet boundary conditions (nothing to be done for linear elements...revise it for quadratic elements ! )'])
[DOFr,dR] = DirichletCONDtime_GENERAL(DIRICHLET,DATA,ndim,MESH,[],OTHER_INPUTS) ;



DISP_CONDITIONS.DOFr = DOFr;
DOFl = 1:DATA.MESH.ndof ;
DOFl(DOFr) = []  ;
DISP_CONDITIONS.DOFl = DOFl ;
DISP_CONDITIONS.dR = dR;
dR = DefaultField(dR,'RIGID_BODY_MOTION',[]) ;
DISP_CONDITIONS.RIGID_BODY_MOTION = dR.RIGID_BODY_MOTION ;
dR.RIGID_BODY_MOTION = [] ;
OPERFE.DOFm = []  ; 

% 7. Nodal body forces due to gravity (SELF-WEIGHT)
% disp(['Nodal body forces due to gravity (self-weight) ---only valid for linear elements, revise it for quadratic '])
%
% gV = repmat(DATA.vGRAVITY(:),nnode,1) ;
% gVall = zeros(DATA.MESH.ndof,1);
% gVall(1:length(gV)) = gV ;
Fbody.U = zeros(DATA.MESH.ndof,1) ;
Fbody.a = ones(size(DATA.STEPS)) ;
OTHER_data.Fbody_U  =Fbody.U ;

% Fbody.a = ones(size(DATA.STEPS)) ;    % Space-Time decomposition.  Gravity is the same for all time steps

% for idim = 1:ndim
%     gV(idim:ndim:end) = gV(idim:ndim:end).*densGLO;
% end
% Fbody.U = NstW'*gV ;
% % %dR = DefaultField(dR,'RIGID_BODY_MOTION',[] ) ;
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
%disp('Imposed forces option not available yet...')
Ftrac.U =  zeros(DATA.MESH.ndof,1);
Ftrac.a = ones(size(DATA.STEPS))  ;

% Ftrac = NeumannCONDtime(NEUMANN,DATA,ndim,MESH,NstT_W_N_boundaries,GEOproperties) ;
% if  ~isempty(Ftrac.U)
%     if any(Ftrac.U)
%         % See explanation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/18_BEAMQ4.mlx
%         % Archetypal file
%         % in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/18_BEAMQ4/INPUTS_EIFE_shapeBUBloc.m
%         Unew = zeros(length(DATA.MESHextended.DOFS_TO_KEEP),size(Ftrac.U,2)) ;
%         Unew(DATA.MESHextended.DOFS_bLOC,:) = Ftrac.U ;
%         Ftrac.U = Unew ;
%     else
%         Ftrac.U = [] ;
%         Ftrac.a = [] ;
%     end
% end

% if ~isempty(DISP_CONDITIONS.RIGID_BODY_MOTION)
%     % Rotation of body forces
%     % -----------------------
%     Ftrac =    RotateForces(Ftrac,DISP_CONDITIONS.RIGID_BODY_MOTION.RotREFP,DATA,MESH) ;
% % end
%
% OTHER_data.NstT_W_N_boundaries = NstT_W_N_boundaries ;

%
%
%
OPERFE.HYDRO = [] ;
%
% if ~isempty(DATA.FOLLOWER_LOADS)
%     switch  DATA.FOLLOWER_LOADS.TYPE
%         case 'HYDROSTATIC'
%
%             DATA.FOLLOWER_LOADS.HYDROSTATIC = DefaultField(DATA.FOLLOWER_LOADS.HYDROSTATIC,'INCLUDE_DYNAMIC_PRESSURE',0) ;
%
%             % Hydrostatic forces
%             [OPERFE.HYDRO,DATA,MESH ]= HYDROST_computeOPER(MESH,DATA)   ;
%             %
%             DATA.FOLLOWER_LOADS.HYDROSTATIC  = DefaultField(DATA.FOLLOWER_LOADS.HYDROSTATIC,'NAME_MESH_WATERLINE',[]) ;
%             if ~isempty(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE)
%                 % WATERLINE MESH
%                 % ------------------
%                 [DATA.WATERLINE_MESH]= ReadMeshFileStr(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE,'READ_MATERIAL_COLUMN',1)  ;
%
%             end
%
%
%         otherwise
%             error('Option not implemented')
%     end
%
% end

DATA = DefaultField(DATA,'EstimateCoefficientsDampingFromNatFrequencies',0) ; 

if DATA.EstimateCoefficientsDampingFromNatFrequencies == 1
    [OPERFE.Ddamp,OTHER_INPUTS.alphaD,OTHER_INPUTS.betaD]= EstimateDdampFromNatFrequencies_EIFEM(OPERFE,DISP_CONDITIONS,DATA) ;
else
    OTHER_INPUTS = DefaultField(OTHER_INPUTS,'alphaD',0) ;
    OTHER_INPUTS = DefaultField(OTHER_INPUTS,'betaD',0) ;
    
    if OTHER_INPUTS.alphaD == 0
        OPERFE.Ddamp = [] ; %
    else
        OPERFE.Ddamp = OTHER_INPUTS.alphaD*OPERHROM.M ;
        if OTHER_INPUTS.betaD ~= 0
            error('Option not implemented')
        end
    end
    
end



% -------------------------
% 9. INITIAL CONDITIONS
% --------------------------

ndof = DATA.MESH.ndof; % prod(size(MESH.COOR)) ;
%DOFl = setdiff(1:ndof,DOFr) ;

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
