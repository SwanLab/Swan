function DATAOUT=    BCeifem_CWeight(DATALOC)
if nargin == 0
    DATALOC = [] ; 
end

 
APPLIED_LOAD =0     ; 
% LENGTH equivalence
% Length one cell = 0.05 m 
%  Effective length 1 cell = 0.0408 m 
% Therefore, we have to multiply the FE result by 0.05/0.0408 = 1.225  
%[XX1,XX2] = ndgrid(disp_x,force) ;
NameParamStudy = mfilename ;

INPUTS_PARAMETERS = [APPLIED_LOAD] ;

% ------------------------%
t0 = 0 ; tEND = 1 ; ntimes = 2; % 50 ;
DATA_STEPS = linspace(t0,tEND,ntimes) ;
% For printing in GID
STEPS_print_FREQ = 1 ; % Only
% Steps to print in GID (SUBSET OF DATA.STEPS)
DATA_PRINT_NSTEPS =  unique([1:STEPS_print_FREQ:length(DATA_STEPS),length(DATA_STEPS)]) ;


% NAME OF THE .MAT FILES
FOLDER = [cd,filesep,'SNAPSHOTS',filesep]  ;
if ~exist(FOLDER)
    mkdir(FOLDER)
end
DATAOUT.FE_VARIABLES_NAMEstore = [FOLDER,filesep,NameParamStudy,'_FEoper','.mat'] ;  

DATAOUT.INPUTS_PARAMETERS = INPUTS_PARAMETERS; 
DATAOUT.t0 = t0; 
DATAOUT.tEND = tEND; 
DATAOUT.DATA_STEPS = DATA_STEPS; 
DATAOUT.DATA_PRINT_NSTEPS = DATA_PRINT_NSTEPS; 
DATAOUT.NameParamStudy = NameParamStudy; 


DATAOUT.SMALL_STRAIN_KINEMATICS= 1;

DATAOUT.DirichletBoundaryConditions = 'DirichBC_CANT' ; 
%DATAOUT.NeumannBoundaryConditions = [] ; 

%DATAOUT.InputDataFile_FE = 'INPUTS_2param' ; 
%DATAOUT.InputDataFile_HROM = 'INPUTS_HROM_2param' ; 
 
DATAOUT.Print_InternalVarStrain=1 ; 
DATAOUT.typePROBLEM='3D' ; 

DATAOUT.vGRAVITY = [0,-9.81,0] ; 

