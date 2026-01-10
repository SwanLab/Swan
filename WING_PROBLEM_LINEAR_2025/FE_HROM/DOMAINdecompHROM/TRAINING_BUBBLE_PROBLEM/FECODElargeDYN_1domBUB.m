function [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] = ...
    FECODElargeDYN_1domBUB(DATAFILE,PARAMETERS)
% FINITE ELEMENT CODE,  adapted to bubble-type (zero work) boundary
% conditions
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/08_IndividualTRAINbub.mlx
% JAHO, 29-Obt-2023 
if nargin ==0
    load('tmp3.mat')
end
TIME_COMPUT=[] ;
CONVERGED = [] ; 
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
OTHER_output.Fbody = Fbody; 
TIME_COMPUT.INPUTS = toc(TIME_COMPUT.INPUTS) ;
disp(['...done in (',num2str(TIME_COMPUT.INPUTS ),' s)'])
disp('*************************************')
DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0) ;
%DATA = DefaultField(DATA,'PRECOMPUTE_ELASTIC_STIFFNESS_MATRIX',0) ; % 9-Feb-2022. 
PARAMETERS = DefaultField(PARAMETERS,'OnlyComputeMassAndStiffnessMatricesSmallStrains',0) ; % 19-Feb-2022. 
OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[]) ; % 25-Feb-2022. 
PARAMETERS = DefaultField(PARAMETERS,'INTERNAL_FORCES_USING_precomputed_Kstiff',0) ; % 24-May-2022. Avoid computing stresses in linear problems 
DATA = DefaultField(DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff',PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff) ; % 24-May-2022. Avoid computing stresses in linear problems 
 



% 
% if (DATA.SMALL_STRAIN_KINEMATICS == 1 && DATA.PRECOMPUTE_ELASTIC_STIFFNESS_MATRIX  == 1) || ...
%         PARAMETERS.OnlyComputeMassAndStiffnessMatricesSmallStrains == 1 || PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff == 1
%     % Determine Stiffness Matrix, small strains 
%     d  = zeros(DATA.MESH.ndof,1) ;    
%     VAR.DISP =zeros(DATA.MESH.ndof,1) ;VARint_n =[] ;    
%     [~,celastST,FgradST,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;        
%     OTHER_output.K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST) ;
%     if PARAMETERS.OnlyComputeMassAndStiffnessMatricesSmallStrains == 1
%         CONVERGED = 0 ; 
%         return 
%     end    
%     if PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff == 1
%       OPERFE.KinternalFORCES_given = OTHER_output.K ; 
%     end
% end




disp('*************************************')
disp('INITIALIZATIONS')
disp('*************************************')
% Variables at time n (and previous time steps if required. Initialization of
% snapshots )
[DATA,VAR,SNAP] = INITIALIZATIONvar_BUBFE(DATA,INICOND)   ;
%

%d = VAR_n.DISP(:,istep-1) ;
TIME_COMPUT.TIME_STEP_LOOP = tic ;
PARAMETERS = DefaultField(PARAMETERS,'ONLY_PRINT_GID',0) ;
if PARAMETERS.ONLY_PRINT_GID == 0
  %  if DATA.ISDYNAMIC == 0
        [DATA,CONVERGED]=NONLINEARstaticLARGE_bubfe(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
   % else
   %     [DATA,CONVERGED]=NONLINEARdynamicLARGE(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
   % end
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
PARAMETERS = DefaultField(PARAMETERS,'PRINT_RIGID_BODY',0) ;
DISP_CONDITIONS = DefaultField(DISP_CONDITIONS,'RIGID_BODY_MOTION',0 ) ;

if PARAMETERS.PRINT_RIGID_BODY == 1
    DATAinpGID.ADDITION_VARIABLE.DISP = DISP_CONDITIONS.RIGID_BODY_MOTION ;
end
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
    
    
    GidPostProcessLARGE_iterations(MESH,DATA,DATAinpGID);
    
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
DATA.PRINT = DefaultField(DATA.PRINT,'SAVE_AS_BINARY_FILE',0)  ;
if DATA.PRINT.SAVE_AS_BINARY_FILE == 1
    GIDres2bin(DATA.PRINT.BASE_FOLDER) ;
end


