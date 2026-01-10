function [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] =  FEcableslargeDYN(DATAFILE,PARAMETERS)
% FINITE ELEMENT CODE, large strains, dynamic regime
% "Cable" element, only tension, 
% JAHO, 22-Jun-2022
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/...
% FIBREGY_PROJECT_2022/04_MOORING_PROBLEMS/01_STATICgrav.mlx
if nargin ==0
    load('tmp.mat')
end
TIME_COMPUT=[] ;
format long g

% READING AND CONSTRUCTING INPUT DATA FOR THE FE CODE
% -----------------------------------------------------------
% Prototypical file, see
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PENDULUM/INPUTS_1.m
disp('*************************************')
disp('Reading and constructing input data')
disp('*************************************')
TIME_COMPUT.INPUTS = tic ;
[MESH,MATPRO,OPERFE,Fbody,...
    Ftrac,DISP_CONDITIONS,INICOND,DATA,OTHER_output] = ...
    feval(DATAFILE,PARAMETERS) ;
OTHER_output.Ftrac = Ftrac; 
TIME_COMPUT.INPUTS = toc(TIME_COMPUT.INPUTS) ;
disp(['...done in (',num2str(TIME_COMPUT.INPUTS ),' s)'])
disp('*************************************')


disp('*************************************')
disp('INITIALIZATIONS')
disp('*************************************')
% Variables at time n (and previous time steps if required. Initialization of
% snapshots )
[DATA,VAR,SNAP] = INITIALIZATIONvarCABLE(DATA,INICOND)   ;
%

%d = VAR_n.DISP(:,istep-1) ;
TIME_COMPUT.TIME_STEP_LOOP = tic ;
PARAMETERS = DefaultField(PARAMETERS,'ONLY_PRINT_GID',0) ;
if PARAMETERS.ONLY_PRINT_GID == 0
    if DATA.ISDYNAMIC == 0
        [DATA,CONVERGED]=NONLINEARstaticCABLES(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
    else
        [DATA,CONVERGED]=NONLINEARdynamicCABLES(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
    end
else
    disp('...Using results already calculated..')
end
TIME_COMPUT.TIME_STEP_LOOP = toc(TIME_COMPUT.TIME_STEP_LOOP) ;
OTHER_output.TIME_COMPUT = TIME_COMPUT;
OTHER_output.DATA = DATA ; % 25-May-2022 ************************
disp(['...done in (',num2str(TIME_COMPUT.TIME_STEP_LOOP ),' s)'])

disp('Postprocess with GID')
% POST-PROCESS WITH GID

DATAinpGID = [] ;
% PARAMETERS = DefaultField(PARAMETERS,'PRINT_RIGID_BODY',0 ) ;
% if PARAMETERS.PRINT_RIGID_BODY == 1
%     DATAinpGID.ADDITION_VARIABLE.DISP = DISP_CONDITIONS.RIGID_BODY_MOTION ;
% end
for icluster = 1:length(DATA.STORE.NSTEPS_CLUSTER)
    if  isempty(DATA.FOLLOWER_LOADS) ||  isempty(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE)
        GidPostProcessLARGE(MESH,DATA,icluster,DATAinpGID);
    else
        GidPostProcessLARGE_waterline(MESH,DATA,icluster,DATAinpGID);
    end
end

if ~isempty(DATA.SNAP_ITER)
    disp('Printing non-converged results')
    disp('___________________________________')
    DATAlocc= [] ;
    
    
    GidPostProcessLARGE_iterations_CABLE(MESH,DATA,DATAinpGID);
    
   % GidPostProcess_Iterations(MESH.COOR,MESH.CN,MESH.TypeElement,DATA.SNAP_ITER.DISP,DATA.MESH.posgp,'NonConverged',DATAlocc);
end



DATA.SNAP_ITER = [] ;
if iscell(DATA.PRINT.NAME_FILE_MSH)
    fff= DATA.PRINT.NAME_FILE_MSH{1} ; 
else
    fff= DATA.PRINT.NAME_FILE_MSH ; 
end
BASE_FOLDER = fileparts(fff) ;
DATA.PRINT = DefaultField(DATA.PRINT,'BASE_FOLDER',BASE_FOLDER)  ;
DATA.PRINT = DefaultField(DATA.PRINT,'SAVE_AS_BINARY_FILE',1)  ;
if DATA.PRINT.SAVE_AS_BINARY_FILE == 1
    GIDres2bin(DATA.PRINT.BASE_FOLDER) ;
end


