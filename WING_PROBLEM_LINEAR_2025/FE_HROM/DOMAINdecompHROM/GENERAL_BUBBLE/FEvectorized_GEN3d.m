function FEvectorized_GEN3d(NAME_INPUT_FILE,TrainingTrajectoryFun,FILE_MATERIAL_DATA,DATALOC)
%%
%
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
% INPUTS
% -------
if nargin == 0
    NAME_INPUT_FILE     ='INPUTS_1elem' ;  'INPUTS_curvedH' ;  'INPUTS_curved' ;
    TrainingTrajectoryFun = 'TrainingTrajectory' ; 
     FILE_MATERIAL_DATA =  'ElasticMaterialDATA1mat';  'ElasticMaterialDATA' ;  
elseif nargin == 1
     TrainingTrajectoryFun = 'TrainingTrajectory' ; 
     FILE_MATERIAL_DATA =   'ElasticMaterialDATA' ;  
elseif nargin == 2
     FILE_MATERIAL_DATA =   'ElasticMaterialDATA' ;  

end


delete('TRAINING.txt')
diary TRAINING.txt
%DATALOC.ONLY_PRINT_GID =0;
RUN_PROBLEM =1;

DATALOC.FILE_MATERIAL_DATA = FILE_MATERIAL_DATA; 

DATAcommon =    feval(TrainingTrajectoryFun,DATALOC) ;   

DATAcommon = DefaultField(DATAcommon,'SMALL_STRAIN_KINEMATICS',1) ; 
%DATAcommon = DefaultField(DATAcommon,'StoreInitialStiffnessMatrix',1) ; 



DATALOC.SMALL_STRAIN_KINEMATICS= DATAcommon.SMALL_STRAIN_KINEMATICS;
 DATALOC.StoreInitialStiffnessMatrix = 1; 

% value of the input parameter.

DATALOC.FileUsedToRunParametricStudy = [cd,filesep,mfilename] ;
DATALOC.DATAcommon = DATAcommon ; 
% LOOP OVER PARAMETRIC SPACE
DATALOC= DefaultField(DATALOC,'IND_PROJECT_LINEAR',[]) ; 
CASES= 1:length(DATAcommon.INPUTS_PARAMETERS) ; 
CASES_nonlinear = setdiff(CASES,DATALOC.IND_PROJECT_LINEAR) ; 

for ILOADSTATES =1:length(DATAcommon.INPUTS_PARAMETERS)
    DATALOC.PARAM  = DATAcommon.INPUTS_PARAMETERS(ILOADSTATES,:) ;
    DATALOC.ILOADSTATES = ILOADSTATES; 
    DATALOC.NameParamStudy =  [DATAcommon.NameParamStudy,'_param_',num2str(ILOADSTATES)]  ;
    if RUN_PROBLEM == 1
        [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] = ...
            FECODElargeDYN(NAME_INPUT_FILE,DATALOC) ;
        if ILOADSTATES == 1
            save(DATA.FE_VARIABLES_NAMEstore,'MESH','MATPRO','OPERFE','DISP_CONDITIONS','OTHER_output') ;
        else
            if ~isempty(DATALOC.IND_PROJECT_LINEAR) && ILOADSTATES == CASES_nonlinear(1)
                save(DATA.FE_VARIABLES_NAMEstore,'MESH','MATPRO','OPERFE','DISP_CONDITIONS','OTHER_output') ;
            end
        end
    else
        % SHOW SVD ERRORS OF EACH SNAPSHOT
        LocalPostProcess(DATALOC.NameParamStudy) ;
    end
end
diary off
