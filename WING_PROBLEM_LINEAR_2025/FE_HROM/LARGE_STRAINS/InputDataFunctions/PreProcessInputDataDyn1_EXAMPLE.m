
% 
% 

% 0) TIME DISCRETIZATION 
% 
t0 = 0 ; tEND = 5 ; ntimes = 10 ; 
DATA.STEPS = linspace(t0,tEND,ntimes) ; 
DATA.STEPS_PRINT = DATA.STEPS(1:end)  ;      

% PARAMETERS RELATED WITH TIME INTREGRATION SCHEME
DATA.ISDYNAMIC = 1 ; 
DATA.INTEGRATION_SCHEME.TYPE ='NEWMARK' ; 
DATA.INTEGRATION_SCHEME.NEWMARK.gamma = 0.5  ; 
DATA.INTEGRATION_SCHEME.NEWMARK.beta = 0.25  ; 

 
% % 1. NAME OF THE MESH AND DATA FILES. COORDINATES, CONNECTIVITIES,  LISTS
% % OF NODES FOR IMPOSING BOUNDARY CONDITIONS
% %--------------------------------------------------------------------------
DATA.NameFileMeshDATA =[cd,filesep,'PREPROCESS/PEND.msh'] ;  % FE mesh of the studied geometry. Remember to load
% problem type: /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/PROBLEMTYPES_GID/PROBLEM_TYPE_SIMPLE.gid % wherein the mesh has been constructed. To generate the files needed by matlab,
 % remember to export both the mesh and the conditions data. (Files >
 % Export > Gid Mesh) and (Files > Export > Calculation File)
% ---------------------------------------------------------------------------
%-------------------------------------------------------------------------------------
% 2. Type of structural problem (plane stress (pstress), plane strain (pstrain), 3D)
% --------------------------------------------------------------
DATA.typePROBLEM = 'pstress' ;
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% 3. Material data  (linear elasticity)
% -----------------------------------------------------------------------------------
imat =1 ; % Index material
% Elasticity matrix 
E = 206.9e9  ; %  MPa, Young's modulus
nu = 0.29; % Poisson's coefficient
% Compliance matrix for an isotropic materials (with all entries, 3D)
% See slides, page 23. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = E/2/(1+nu) ;  % Shear modulus
celasINV3D = [1/E  -nu/E -nu/E  0 0 0  
           -nu/E  1/E  -nu/E  0 0 0
           -nu/E  -nu/E  1/E  0 0 0           
            0      0    0  1/G   0 0
            0      0      0  0 1/G 0  
            0      0      0  0  0  1/G] ; 
ElasticityMatrix = inv(celasINV3D) ;  
PROPMAT(imat).ElasticityMatrix =  ElasticityMatrix  ; %   
PROPMAT(imat).Density = 7850 ; % Kg/m3

% imat = 2 .... Repeat the above sequence of operations for defining
% another material 

% % -----------------------------------------------------------
% % 4. Dirichlet boundary conditions (prescribed displacements)
% % -----------------------------------------------------------

icond = 1; % Number of condition (fixed along time )
%---------------------------------------------------------------
DIRICHLET(icond).NUMBER_SURFACE = 1 ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
iloadstate = 1;  % Loading stage 
DIRICHLET(icond).PRESCRIBED_DISP = [] ;
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE = {[],[]} ;  % Constraints x,y and z directions. If empty, no constraints  
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t ;  % Function definining the temporal evolution of the constraint
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [0,tEND/2];  % Outside this interval, the amplitude is assumed to be zero

iloadstate = 2; 
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE = {[],[]} ;  %  
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t ;  % 
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [tEND/2,tEND];  %  

 
icond = 2; % Number of condition (fixed along time )
%---------------------------------------------------------------
DIRICHLET(icond).NUMBER_POINT_MESH = 314 ;   % Number of node in which DISPLACEMENT  is prescribed
iloadstate = 1; 
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE = {0,0} ;  %  
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t ;  % 
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [0,tEND] ;  %  


% % -------------------------------------------------
% % 5. External forces  (NEUMANN CONDITIONS, POINT FORCES AND BODY FORCES)
% --------------------------------------------------------------------
%  5.1) NEUMANN COND. Loads per unit surface  
% ------------------------------------------------------------------ 
icond= 1 ;
NEUMANN(icond).NUMBER_SURFACE = 1 ;  % Surface on which the load is applied 
iloadstate = 1; 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE = [] ; 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).AMPLITUDE =  [0,0] ; % Force per unit surface 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).TIMEFUN =  @(t) t ; ; % 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).INTERVAL = [0,tEND/2] ; ; % 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).ISLOCAL = 1 ; ; %  locaL -- [NORMAL COMP., TANGENT. COMP. ]


iloadstate = 2; 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).AMPLITUDE =  [0,0] ; %  
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).TIMEFUN =  @(t) t ; ; %  
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).INTERVAL = [tEND/2,tEND] ; ;  
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).ISLOCAL = 1 ; %  
 
icond= 2;
NEUMANN(icond).NUMBER_SURFACE = 2 ;  % Surface on which the load is applied
iloadstate = 1; 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).AMPLITUDE =  [0,0] ; %  
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).TIMEFUN =  @(t) t ; ; % 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).INTERVAL = [0,tEND] ; ; % 
NEUMANN(icond).FORCE_PER_UNIT_SURFACE(iloadstate).ISLOCAL = 1 ; ; %   

icond= 3;
NEUMANN(icond).NUMBER_NODE = [2  3 4];  % Number of finite element node in which the force is applied 
iloadstate = 1; 
NEUMANN(icond).FORCE(iloadstate).AMPLITUDE =  [0,0] ; %  
NEUMANN(icond).FORCE(iloadstate).TIMEFUN =  @(t) t ; ; %  
NEUMANN(icond).FORCE(iloadstate).INTERVAL = [0,tEND] ; ; % 
% 
 
% % --------------------------------------------------------------------
% %---5.3)  Body forces ***** only  GRAVITY ****
% % ----------------------------------------------------------------------
DATA.vGRAVITY =[0,-9.81];  % Gravity vector 

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
% 
% [MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATA] ...
%     =  PreProcessInputDataDyn1(DATA,PROPMAT,DIRICHLET,NEUMANN,INITIAL_CONDITIONS) ; 
%  