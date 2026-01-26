function EIFEvectorizedGEN2(NAME_INPUT_FILE,TrainingTrajectoryFun,FILE_MATERIAL_DATA,DATALOC)
%%
%% ---------------------------------------------------------------------------------------
if ~exist('FECODElargeDYN')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
% INPUTS
% ----------------------------------------------------------------------------------------
% -------
if nargin == 0
    NAME_INPUT_FILE     ='INPUTS_1elem' ; 
    TrainingTrajectoryFun = 'TrainingTrajectory' ;
    FILE_MATERIAL_DATA =  'ElasticMaterialDATA1mat';  
elseif nargin == 1
    TrainingTrajectoryFun = 'TrainingTrajectory' ;
    FILE_MATERIAL_DATA =   'ElasticMaterialDATA' ;
elseif nargin == 2
    FILE_MATERIAL_DATA =   'ElasticMaterialDATA' ;
    DATALOC = [] ; 
elseif nargin == 3
    DATALOC = [] ; 
end


delete('TRAININGeife.txt')
diary TRAININGeife.txt
DATALOC = DefaultField(DATALOC,'ONLY_PRINT_GID',0) ;% =0;
RUN_PROBLEM =1;

DATALOC.FILE_MATERIAL_DATA = FILE_MATERIAL_DATA;

DATAcommon =    feval(TrainingTrajectoryFun,DATALOC) ;
DATALOC.SMALL_STRAIN_KINEMATICS= 1;


% value of the input parameter.

DATALOC.FileUsedToRunParametricStudy = [cd,filesep,mfilename] ;
DATALOC.DATAcommon = DATAcommon ;
% LOOP OVER PARAMETRIC SPACE
for ILOADSTATES =1:size(DATAcommon.INPUTS_PARAMETERS(:),2)
    DATALOC.PARAM  = DATAcommon.INPUTS_PARAMETERS(ILOADSTATES,:) ;
    
    DATALOC.NameParamStudy =  [DATAcommon.NameParamStudy,'_param_',num2str(ILOADSTATES)]  ;
    if RUN_PROBLEM == 1
        [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,DISP_CONDITIONS,OTHER_output] = ...
            EIFECODElargeDYN(NAME_INPUT_FILE,DATALOC) ;
        if ILOADSTATES == 1 && ~isempty(MATPRO)
            save(DATA.FE_VARIABLES_NAMEstore,'MESH','MATPRO','OPERFE','DISP_CONDITIONS','OTHER_output') ;
        end
    else
        % SHOW SVD ERRORS OF EACH SNAPSHOT
        LocalPostProcess(DATALOC.NameParamStudy) ;
    end
end
diary off
