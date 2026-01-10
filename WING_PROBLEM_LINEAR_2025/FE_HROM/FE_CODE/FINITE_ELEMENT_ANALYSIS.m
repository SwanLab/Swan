function  FINITE_ELEMENT_ANALYSIS(NAME_DATA_INP,DATAINP)
% Linear/Nonlinear FE analysis
% ........................
% INPUTS:  NAME_DATA = Name input data file 
%          DATAINP.EXECUTABLE_FOLDER = Folder in which the FE routires are located 
%          DATAINP.FOLDER_INPUT_FILE = Folder  in which NAME_DATA is
%          located
% EXAMPLE of Input data files:  '/home/joaquin/Desktop/CURRENT_TASKS/DOCTORAL_THESIS_SUPERVISED/AGUSTINA_ANACONDA/DOCS/SIMULATIONS/ElastCellSil/'
% 
%  
%  
% 22-Jun-2019. JAHO
% -----------------------------------------------------------------------------

if nargin == 0
    NAME_DATA_INP ='CylinderPressurePSTRAIN' ; 'CylinderPressurePSTRAINthin' ;
    DATAINP.FOLDER_INPUT_FILE = '/home/joaquin/Desktop/CURRENT_TASKS/DOCTORAL_THESIS_SUPERVISED/AGUSTINA_ANACONDA/DOCS/SIMULATIONS/ElastCellSil'; 
elseif  nargin == 1
        DATAINP.FOLDER_INPUT_FILE = '/home/joaquin/Desktop/CURRENT_TASKS/DOCTORAL_THESIS_SUPERVISED/AGUSTINA_ANACONDA/DOCS/SIMULATIONS/'; 

end

CURRENT_FOLDER = cd ; 

EXECUTABLE_FOLDER_def = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/' ;  
DATAINP = DefaultField(DATAINP,'EXECUTABLE_FOLDER',EXECUTABLE_FOLDER_def) ; 

EXECUTABLE_FOLDER = DATAINP.EXECUTABLE_FOLDER ; 

% Adding folders
try 

NAME_DATA = [ DATAINP.FOLDER_INPUT_FILE,filesep,NAME_DATA_INP] ; 

run([NAME_DATA]) ;
catch 
    DATAINP.FOLDER_INPUT_FILE = CURRENT_FOLDER ; 
    NAME_DATA = [ DATAINP.FOLDER_INPUT_FILE,filesep,NAME_DATA_INP] ; 

    run([NAME_DATA]) ;

end
if ~exist('GeometryRVE','file') ; addpath('FE_CODE/MULTI_ROM/'); end
if ~exist('FE_ELASTOSTATIC','file') ; addpath([EXECUTABLE_FOLDER,'FE_CODE']) ; end
if ~exist('FE_NONLINEAR','file') ; addpath([EXECUTABLE_FOLDER,'FE_CODE',filesep,'NONLINEAR_Dyn']) ; end

% ----------------------------------
% --------------------------------------------------------
% 3. MATERIAL PROPERTIES  (elastic)
% ----------------------
FUNinput.INPUTS.MATERIAL =  NAME_INPUT_DATA_MATERIAL   ;
% 4. Other inputs
cd(EXECUTABLE_FOLDER) ;
run(NAME_INPUT_DATA_MATERIAL) ;
cd(EXECUTABLE_FOLDER) ;

% Definition of traction forces
% ----------------------------

FUNinput.NAME = 'INPUTS_GENERAL_NONLINEAR' ; % NAME OF THE FUNCTION THAT PREPARES THE INPUT DATA FOR THE SOLVER
% File that constructs the vector of nodal tracion forces:
FUNinput.INPUTS.NEUMANN_BOUNDARY_CONDITIONS = 'NEUMANN_CONDITIONS_GID_PROBLEMTYPE' ;
if ~exist('FORCES','var')
    FORCES = []  ;
end
FUNinput.INPUTS.FORCES  = FORCES;


DATA.INPUTDATAfile = [FOLDER,filesep,NAME_DATA] ;
% Dirichlet Boundary Conditions
FUNinput.INPUTS.DISPLACEMENT_COND_FILE ='DIRICHLET_CONDITIONS_GID_PROBLEMTYPE' ;
FUNinput.INPUTS.DISP  = DISP;



FUNinput.INPUTS.NameFileMesh = [FOLDER,filesep,NAME_MESH] ; ;
%   FUNinput.INPUTS.NAMEPROJECT = NAMES_LOC_PROJ{iprojects} ;



DATA.RECALCULATE_STIFFNESS =1 ;  % I disable this to avoid conflicts (27-Sept-2018)%
% Other variables
DATA.CALCULATE_MASSMATRIX =0 ;  % Compute mass matrix
DATA.STORE_STIFFNESS = 0; %Store stiffness matrix and other variables
DATA_TYPESOLVER = 0; % Conjugated gradient if = 1.
DATA_niterCONJG = 100000 ; % Number of iterations  conjugated gradient
% DEFAULT INPUTS ---------
% -----------------------------------------


AAAA = tic ;

DATA = DefaultField(DATA,'TIME_DISCRETIZATION',[]) ; 
if isempty(DATA.TIME_DISCRETIZATION)
    
    % Linear problems 
    FUNinput.NAME = 'INPUTS_GENERAL' ; % NAME OF THE FUNCTION THAT PREPARES THE INPUT DATA FOR THE SOLVER

    DATA.StrainStressWith4Components = 1 ; 
    FE_ELASTOSTATIC(FUNinput,DATA) ;
else
    % Nonlinear problem. J2 plasticity
    FE_NONLINEAR(FUNinput,DATA) ;
end


AAAA = toc(AAAA) ;
disp(['Total time FE training simulation = ',num2str(AAAA),' s'])

cd(FOLDER) ;
