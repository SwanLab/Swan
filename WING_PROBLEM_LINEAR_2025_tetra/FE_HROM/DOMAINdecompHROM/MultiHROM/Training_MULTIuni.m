function Training_MULTIuni(INPUT_PROBLEMS,DATAoffline) ; 
%% TRAINING FUNCTION, MULTISCALE APPROACH 
% SEE COMMENTS IN /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
if nargin == 0
    INPUT_PROBLEMS ='DefineInputs_HOMOG' ; % 'DefineInputs_3DLIN'; 'DefineInputs_HOMOG_quad' ;   % FUNCTION DEFINING  SPECIFIC INPUTS OF THIS PROBLEM  
elseif nargin == 1
    DATAoffline = [] ; 
    
end
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
% INPUTS
% -------
delete('TRAINING.txt') ; diary TRAINING.txt ; 
%-------------------------------------------------
DATALOC =    feval(INPUT_PROBLEMS) ;    % Recovering the variables defined inside function "INPUT_PROBLEMS" 
DATAoffline = DefaultField(DATAoffline,'NAMEMESH', DATALOC.NameFileMeshDATA); 
DATALOC.NameFileMeshDATA = DATAoffline.NAMEMESH ; 

%DATAcommon =    feval(TrainingTrajectoryFun) ;   


DATALOC.FileUsedToRunParametricStudy = [cd,filesep,mfilename] ; 
% LOOP OVER PARAMETRIC SPACE
Fbody = cell(1,size(DATALOC.INPUTS_PARAMETERS,2) ) ; 
Ftrac = cell(1,size(DATALOC.INPUTS_PARAMETERS,2) ) ; 

for ILOADSTATES =1:length(DATALOC.INPUTS_PARAMETERS) 
    DATALOC.PARAM  = DATALOC.INPUTS_PARAMETERS{ILOADSTATES} ;    
    DATALOC.NameParamStudy =  [DATALOC.NameParamStudyLOC,'_param_',num2str(ILOADSTATES)]  ;
        [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] = ...
            FECODElargeDYN(DATALOC.InputDataFile_FE ,DATALOC) ;
        if ILOADSTATES == 1
            save(DATA.FE_VARIABLES_NAMEstore,'MESH','MATPRO','OPERFE','DISP_CONDITIONS','OTHER_output') ;
            NAME_FE_WS_1st = DATA.FE_VARIABLES_NAMEstore ; 
        end
        Fbody{ILOADSTATES} = OTHER_output.Fbody ; 
        Ftrac{ILOADSTATES} = OTHER_output.Ftrac ;
end
save(NAME_FE_WS_1st,'Fbody','Ftrac','-append') ;
diary off
