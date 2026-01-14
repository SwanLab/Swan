function Training_MULTIfluct(INPUT_PROBLEMS,DATAoffline) ; 
%% TRAINING FUNCTION, MULTISCALE APPROACH 
% SEE COMMENTS IN /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% New version, store fluctuations 
% JAHO, 3-Dec-2023
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/14_GIVEN_FLUCT_Q4.mlx
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
dBoundary_FLUCT =  cell(1,size(DATALOC.INPUTS_PARAMETERS,2) ) ; % Fluctuations 

for ILOADSTATES =1:size(DATALOC.INPUTS_PARAMETERS,2) 
    DATALOC.PARAM  = DATALOC.INPUTS_PARAMETERS(:,ILOADSTATES) ;    
    DATALOC.NameParamStudy =  [DATALOC.NameParamStudyLOC,'_param_',num2str(ILOADSTATES)]  ;
        [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] = ...
            FECODElargeDYN(DATALOC.InputDataFile_FE ,DATALOC) ;
        if ILOADSTATES == 1
            save(DATA.FE_VARIABLES_NAMEstore,'MESH','MATPRO','OPERFE','DISP_CONDITIONS','OTHER_output') ;
            NAME_FE_WS_1st = DATA.FE_VARIABLES_NAMEstore ; 
        end
        Fbody{ILOADSTATES} = OTHER_output.Fbody ; 
        Ftrac{ILOADSTATES} = OTHER_output.Ftrac ;
        % Extract fluctuations 
        % -----------------------
       dBoundary_FLUCT{ILOADSTATES}=  ExtractFluctuationsSolidEIFEM(INPUT_PROBLEMS,DATAoffline,ILOADSTATES,...
            MESH,MATPRO,OPERFE,DISP_CONDITIONS,OTHER_output) ; 
        
end

 Ufluct =  PrepareFluctuationsEIFEM(INPUT_PROBLEMS,DATAoffline,dBoundary_FLUCT,...
            MESH,MATPRO,OPERFE,DISP_CONDITIONS,OTHER_output) ; 
        
        
        




save(NAME_FE_WS_1st,'Fbody','Ftrac','-append') ;

% TRAINING WITH NEW BOUNDARY CONDITIONS (PRESCRIBED FLUCTUATIONS)
% -------------------------------------------
DATALOC.Ufluct = Ufluct ; 
for ILOADSTATES =1:size(DATALOC.INPUTS_PARAMETERS,2) 
    DATALOC.PARAM  = DATALOC.INPUTS_PARAMETERS(:,ILOADSTATES) ;    
    DATALOC.NameParamStudy =  [DATALOC.NameParamStudyLOC,'_param_',num2str(ILOADSTATES)]  ;
        [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] = ...
            FECODElargeDYN(DATALOC.InputDataFile_FE ,DATALOC) ;
        if ILOADSTATES == 1
            save(DATA.FE_VARIABLES_NAMEstore,'MESH','MATPRO','OPERFE','DISP_CONDITIONS','OTHER_output') ;
            NAME_FE_WS_1st = DATA.FE_VARIABLES_NAMEstore ; 
        end
        Fbody{ILOADSTATES} = OTHER_output.Fbody ; 
        Ftrac{ILOADSTATES} = OTHER_output.Ftrac ;
        % Extract fluctuations 
        % -----------------------
       dBoundary_FLUCT{ILOADSTATES}=  ExtractFluctuationsSolidEIFEM(INPUT_PROBLEMS,DATAoffline,ILOADSTATES,...
            MESH,MATPRO,OPERFE,DISP_CONDITIONS,OTHER_output) ; 
        
end

save(NAME_FE_WS_1st,'Fbody','Ftrac','-append') ;



diary off
