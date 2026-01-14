function [MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,...
    INICOND,DATA,OTHER_output] =  INPUTS_FE_LOC(DATAINPUT)
% % See README_4paramINP.pdf,
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DOCS/06_ReducedOrderModel_01.pdf
% THIS FUNCTION PREPARES THE DATA NEEDED TO RUN THE FE LARGE STRAIN
%   FUNCTION .....PreProcessInputDataDyn1... 
%  
% ANALYSIS.  
% INPUTS
%
% PARAMETERS

% OUTPUTS
% 1) MESH --- Array structure with coordinates, connectivities....  sdsd
% 2) MATPRO: Matrices and other variables characterizing material properties (defined at integration points)
%            for instance, the elasticity matrix
%
% 3) OPERFE: Miscellaneous operators used directly in the formulation (mass matrix, for instance)
% 4) Fbody.  Body forces  (in a space-time separated fashion)
% 5) Ftrac.  Forces applied on the surface  (in a space-time separated fashion)
% 6) DISP_CONDITIONS: Displacement boundary conditions
% 7) INICOND: Initial conditions (displacement and velocity)
% 8) DATA. Miscellaneous data (integration scheme, name of the mesh...)

% 0) TIME DISCRETIZATION
% % ------------------------%
 t0 = DATAINPUT.DATAcommon.t0 ; tEND =  DATAINPUT.DATAcommon.tEND ;  
 DATA.STEPS = DATAINPUT.DATAcommon.DATA_STEPS ;
% % For printing in GID
% STEPS_print_FREQ = 1 ; % Only
% % Steps to print in GID (SUBSET OF DATA.STEPS)
 DATA.PRINT.NSTEPS =   DATAINPUT.DATAcommon.DATA_PRINT_NSTEPS ; 


% PARAMETERS RELATED WITH TIME INTREGRATION SCHEME
DATA.ISDYNAMIC = 0 ;
DATA.INTEGRATION_SCHEME.TYPE ='NEWMARK' ;
DATA.INTEGRATION_SCHEME.NEWMARK.gamma = 0.5  ;
DATA.INTEGRATION_SCHEME.NEWMARK.beta = 0.25  ;
DATA.INTEGRATION_SCHEME.steps_KINEM_VARIABLES_TO_STORE = 1; % Number of previous time steps information
% to be stored  (with a view towards extending to other time-integrating systems
% (or for doing extrapolations using, say, the SVD...). By default = 1 )% figure(1)
% hold on
 

DATA.INTEGRATION_SCHEME.steps_STRESS_VARIABLES_TO_STORE =1; 


DATA.SMALL_STRAIN_KINEMATICS = DATAINPUT.SMALL_STRAIN_KINEMATICS ; 



% TOLERANCE NEWTON_RAPHSON
DATA.NEWTON_RAPHSON.TOL_FORCES_REL = 1e-6 ; %1e-8;  
DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS = 1e-10;  

DATA.NEWTON_RAPHSON.NMAXiter = 50;  


% % 1. NAME OF THE MESH AND DATA FILES. COORDINATES, CONNECTIVITIES,  LISTS
% % OF NODES FOR IMPOSING BOUNDARY CONDITIONS
% %--------------------------------------------------------------------------
NAME ='FEpore' ;  'beam_2x2_rep16_fineINCq' ; 
NameFileMeshDATA =[cd,filesep,'GIDPRE/',NAME,'.msh'] ; 

DATAINPUT  = DefaultField(DATAINPUT,'NameFileMeshDATA',NameFileMeshDATA); 

DATA.NameFileMeshDATA = DATAINPUT.NameFileMeshDATA ; 


%[cd,filesep,'GIDPRE/RECT_100_quad.msh'] ;    % FE mesh of the studied geometry. Remember to load
% problem type: /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/PROBLEMTYPES_GID/PROBLEM_TYPE_SIMPLE.gid % wherein the mesh has been constructed.
% To generate the files needed by matlab,
% remember to export both the mesh and the conditions data. (Files >
% Export > Gid Mesh) and (Files > Export > Calculation File)
% ---------------------------------------------------------------------------
%-------------------------------------------------------------------------------------
% 2. Type of structural problem (plane stress (pstress), plane strain (pstrain), 3D)
% --------------------------------------------------------------
DATAINPUT.DATAcommon  = DefaultField(DATAINPUT.DATAcommon,'typePROBLEM','pstrain'); 

DATA.typePROBLEM = DATAINPUT.DATAcommon.typePROBLEM ;
% -----------------------------------------------------------------------------------

%[PROPMAT,DATA ]= ElasticMaterialDATA(DATA) ; 
DATAINPUT = DefaultField(DATAINPUT,'FILE_MATERIAL_DATA','ElasticMaDATA2mat') ; 
[PROPMAT,DATA ]= feval(DATAINPUT.FILE_MATERIAL_DATA,DATA);
%DATALOC.SMALL_STRAIN_KINEMATICS= 1;
 


% BOUNDARY CONDITIONS FOR DISPLACEMENTs
[DIRICHLET] = feval(DATAINPUT.DATAcommon.DirichletBoundaryConditions,DATAINPUT) ;  % DirichletBoundaryConditionsFUN(tEND) ; 
% figure(1)
% hold on
% xlabel('Time')
% ylabel('Amp disp')
% for iloadstate = 1:length(DIRICHLET(icond).PRESCRIBED_DISP)
%     fplot(DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN,DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL)
% end


 

% % -------------------------------------------------
% % 5. External forces  (NEUMANN CONDITIONS, POINT FORCES AND BODY FORCES)
if isfield(DATAINPUT.DATAcommon,'NeumannBoundaryConditions')  
[NEUMANN] = feval(DATAINPUT.DATAcommon.NeumannBoundaryConditions,DATAINPUT) ;  % DirichletBoundaryConditionsFUN(tEND) ;
else
NEUMANN = [] ; 
end
% %

% % --------------------------------------------------------------------
% %---5.3)  Body forces ***** only  GRAVITY ****
% % ----------------------------------------------------------------------


switch DATA.typePROBLEM
    case '3D'
        vGRAVITY =0*[0,0,-9.81];  % Gravity vector
    otherwise
        vGRAVITY =0*[0,-9.81];  % Gravity vector
        
end

DATAINPUT.DATAcommon = DefaultField(DATAINPUT.DATAcommon,'vGRAVITY',vGRAVITY) ; 

DATA.vGRAVITY = DATAINPUT.DATAcommon.vGRAVITY ; 


% 6) INITIAL CONDITIONS
% --------------------------
% DISPLACEMENTS
%--------------------------
INITIAL_CONDITIONS.DISP.ISZERO =  1;
%INITIAL_CONDITIONS.DISP.ISZERO =  0;   ---> Then, load it from a .mat
%file

% VELOCITY
% -------------------------
INITIAL_CONDITIONS.VELOC.ISZERO =  1;
%INITIAL_CONDITIONS.VELOC.ISZERO =  0;   ---> Then, load it from a .mat
%FILE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% END INPUT DATA


% PreProcessInputDataDyn1 is defined in
% / MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions

OTHER_INPUTS = [] ;
DATAINPUT.DATAcommon = DefaultField(DATAINPUT.DATAcommon,'TypeFunctionDisplacementTRAINING',''); 

OTHER_INPUTS.TypeFunctionDisplacementTRAINING = DATAINPUT.DATAcommon.TypeFunctionDisplacementTRAINING ; 

%
[MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATA,OTHER_output] ...
    =  PreProcessInputDataDyn1(DATA,PROPMAT,DIRICHLET,NEUMANN,INITIAL_CONDITIONS,OTHER_INPUTS) ;


DATA.MESH.nnode = size(MESH.COOR,1) ;
DATA.MESH.ndim = size(MESH.COOR,2) ;
DATA.MESH.nelem = size(MESH.CN,1) ;

DATA.MESH.ndof = DATA.MESH.nnode*DATA.MESH.ndim ;
if ~isfield(DATA.MESH,'nstrain')
    DATA.MESH.nstrain =size(MATPRO.celasglo,2);
end


% ---------------------------------
% STORING INFORMATION
% ----------------------------------
%%%% DECIDING HOW MANY SNAPSHOTS TO STORE
% Notice that we can compute at this level memory requirements to store
% A)  displacements  Snapshots
% B)  Strain, stress snapshots....
% The idea is to avoid storing in the RAM, but rather group some snapshots
% and write it them to disk (after applying the SVD ). Likewise, as this is
% done, it is possible to determine a basis matrix for the corresponding
% variable.

% Assume we have 5 Gb available. Matrices over, say, 100 Mb should be
% written to disk. This is the limit
LIMIT = 50 ; % ;Mb
% What is the size of nodal variables
%
nsizeNOD = DATA.MESH.ndof*length(DATA.STEPS)*8*1e-6 ;
nsizeGAUSS = DATA.MESH.ndofSTRESS*length(DATA.STEPS)*8*1e-6 ;

nclusters = max(ceil(nsizeNOD/LIMIT),ceil(nsizeGAUSS/LIMIT)) ;

%nclusters = 2; 


DATA.STORE.NSTEPS_CLUSTER = cell(1,nclusters) ;
% We have to divide  1:length(DATA.STEPS) into nclusters
FREQ  = ceil(length(DATA.STEPS)/nclusters) ;
iini = 1;
%ifin = FREQ ;
for i=1:nclusters
    ifin = iini + FREQ-1 ;
    ifin = min(ifin,length(DATA.STEPS)) ;
    DATA.STORE.NSTEPS_CLUSTER{i} = iini:ifin ;
    iini = ifin +1 ;
end

DATA.STORE.NAME_MATFILE_STORE = cell(1,nclusters);
DATA.STORE.VAR.VEL= 0;
DATA.STORE.VAR.DISP = 1;
DATA.STORE.VAR.PK2STRESS = 0;
DATA.STORE.VAR.GLSTRAINS = 0;
DATA.STORE.VAR.CAUCHY_STRESS = 0;
DATA.STORE.VAR.RESID = 1;
DATA.STORE.VAR.FEXT = 1;
 DATA.STORE.VAR.STRAIN_ENERGY = 1;

 

 %DATA.STORE.VAR.PK2STRESS =cell(1,nclusters);
% DATA.STORE.VAR.GLSTRAINS =cell(1,nclusters);

% NAME OF THE .MAT FILES
FOLDER = [cd,filesep,'SNAPSHOTS',filesep]  ;
if ~exist(FOLDER)
    mkdir(FOLDER)
end

% SAVE MESH, MATERIAL AND FE OPERATORS 
% ------------------------------------

DATA.FE_VARIABLES_NAMEstore = [FOLDER,filesep,DATAINPUT.NameParamStudy,'_FEoper','.mat'] ;  


MATBASE_FOLDER = [FOLDER,filesep,DATAINPUT.NameParamStudy,filesep] ;

if exist(MATBASE_FOLDER)
    rmdir(MATBASE_FOLDER,'s') ;
end
mkdir(MATBASE_FOLDER)

INFO_SNAPSHOTS.FileUsedToRunParametricStudy = DATAINPUT.FileUsedToRunParametricStudy ;
INFO_SNAPSHOTS.InputFileFunction = [cd,filesep,mfilename] ;
INFO_SNAPSHOTS.InputFileFunction_data = DATAINPUT;

INFO_SNAPSHOTS.EMPLOYED_INPUTS.PROPMAT = PROPMAT;
INFO_SNAPSHOTS.EMPLOYED_INPUTS.DIRICHLET = DIRICHLET;
INFO_SNAPSHOTS.EMPLOYED_INPUTS.NEUMANN = NEUMANN;
INFO_SNAPSHOTS.EMPLOYED_INPUTS.INITIAL_CONDITIONS = INITIAL_CONDITIONS;

INFO_SNAPSHOTS.DATE = datetime ; % Data and time. To check the version employed in the program

name_BASIC_INFO = [MATBASE_FOLDER,'INFO_SNAPSHOTS','.mat'] ;

%FFF =fieldnames(DATA.STORE.VAR);

%for iii = 1:length(FFF)
 %   mkdir([MATBASE_FOLDER,filesep,FFF{iii},filesep])  ;
  %  MATBASE_FOLDER_VAR = [MATBASE_FOLDER,filesep,FFF{iii},filesep] ;
    for icluster = 1:nclusters
        DATA.STORE.NAME_MATFILE_STORE{icluster} = [MATBASE_FOLDER,'SNAP_',num2str(icluster),'.mat'] ;
    end
%end

% COMPRESSING THE VARIABLES...
DATA.STORE.COMPRESS_WITH_SVD = 1; %% SNAPSHOTS ARE DECOMPOSED AS X = U*S*V'
DATA.STORE.TOLERANCE_SVD_COMPRESSION = 1e-6; % Compression tolerance  for the SVD

% ------------------
% PRINTING IN GID
% -----------------

% NAME OF THE .MAT FILES
FOLDER = [cd,filesep,'GIDPOST',filesep]  ;
if ~exist(FOLDER)
    mkdir(FOLDER)
end

NAMERES_folder=  [FOLDER,filesep,DATAINPUT.NameParamStudy,filesep] ;

if ~exist(NAMERES_folder)
    mkdir(NAMERES_folder)
end
NAMERES = [NAMERES_folder,DATAINPUT.NameParamStudy] ; 

DATA.PRINT.NAME_FILE_MSH   =  cell(1,nclusters) ;
DATA.PRINT.NAME_FILE_RES   =  cell(1,nclusters) ;
for icluster = 1:nclusters
    nameloc = [NAMERES,'_clus_',num2str(icluster)] ;
    DATA.PRINT.NAME_FILE_MSH{icluster} = [nameloc,'.msh'] ;
    DATA.PRINT.NAME_FILE_RES{icluster} = [nameloc,'.res'] ;
end

%-------------------------
% Variables to print
% ------------------------
% INFORMATION FOR PRINTING VARIABLES IN GID
% ----------------------------------
NODESV_PROP.DISP.PRINT= 1;
NODESV_PROP.ACEL.PRINT =  0 ;
NODESV_PROP.VEL.PRINT =  0 ;
NODESV_PROP.RESID.PRINT = 0;
NODESV_PROP.FEXT.PRINT = 1;


% INFORMATION FOR PRINTING VARIABLES BY GID
% -----------------------------------------
GAUSSV_PROP.PK2STRESS.PRINT = 0 ;
GAUSSV_PROP.GLSTRAINS.PRINT = 0 ;


GAUSSV_PROP.VONMISES_CAUCHY_STRESS.PRINT = 1 ;
DATA.STORE.VAR.VONMISES_CAUCHY_STRESS = 1;
DATA.PRINT.SAVE_AS_BINARY_FILE  = 0; 

switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
    case 'SMALL_STRAINS_J2_PLASTICITY'
        if DATAINPUT.DATAcommon.Print_InternalVarStrain ==1 
        GAUSSV_PROP.InternalVarStrain.PRINT = 1 ;
        DATA.STORE.VAR.InternalVarStrain = 1;
        end
end


GAUSSV_PROP.CAUCHY_STRESS.PRINT = 0 ;



% ---------------------------------------------------------------------
DATA.PRINT.GAUSSV_PROP = GAUSSV_PROP ;
DATA.PRINT.NODESV_PROP = NODESV_PROP ;

%%%%%%%%%%%%%%%%%%%%%%

INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA = DATA;

save(name_BASIC_INFO,'INFO_SNAPSHOTS')  ;












%
% switch DATA.TypeImplementation
%     case 'J2_plasticity_small_strains'
%         GAUSSV_n.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
%         GAUSSV_n.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
%         GAUSSV_n.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
%
%         GAUSSV_PROP.alpha.TYPE = 'Scalar' ;
%         GAUSSV_PROP.alpha.PRINT = 1 ;  % GID's printing
%         GAUSSV_PROP.alpha.LEGEND = 'INT.VAR.' ;
%
%         GAUSSV_n.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
%         GAUSSV_n.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
%         GAUSSV_n.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
%         GAUSSV_n.VonMises = zeros(OPERfe.ngaus,1);; % stressST-like internal variable
%
%
%
%
%         GAUSSV_PROP.alpha.TYPE = 'Scalar' ;
%
%         DATA.PRINT_GID = DefaultField(DATA.PRINT_GID,'INTERNAL_VARIABLES',1) ;
%         if DATA.PRINT_GID.INTERNAL_VARIABLES == 1
%             GAUSSV_PROP.alpha.PRINT = 1 ;  % GID's printing
%         else
%             GAUSSV_PROP.alpha.PRINT = 0 ;
%         end
%         GAUSSV_PROP.alpha.LEGEND = 'INT.VAR.' ;
%
%
%         DATA.STORE = DefaultField(DATA.STORE,'INTERNAL_VARIABLES',1) ;
%         if DATA.STORE.INTERNAL_VARIABLES == 0
%             GAUSSV_PROP.alpha.STORE = 0 ;
%         end
%
%
%         GAUSSV_PROP.VonMises.TYPE = 'Scalar' ;
%         GAUSSV_PROP.VonMises.PRINT = 1 ;  % GID's printing
%         GAUSSV_PROP.VonMises.LEGEND = 'Von Mises' ;
%     otherwise
%         error('Option not implemented')
% end
%
%
%
% % -----------------------------------------------------------------------
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GAUSSV_np1.EP   = zeros(OPERfe.ngausT,1); % Plastic deformation at time n
% % GAUSSV_np1.stressST = zeros(OPERfe.ngausT,1);  % stressSTes at time n+1
% % GAUSSV_np1.alpha   = zeros(OPERfe.ngaus,1);  % Strain-like internal variable at time n
% % GAUSSV_np1.sigmay = PROPMAT.sigmay_0; % stressST-like internal variable
% % % -----------------------------------------------------------------------



% All the input information is stored in the folder so that, eventually,
% one can recover the exact data used to generate the snapshots

%