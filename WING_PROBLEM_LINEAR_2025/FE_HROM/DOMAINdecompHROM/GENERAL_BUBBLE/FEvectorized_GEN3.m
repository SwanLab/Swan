function FEvectorized_GEN3(NAME_INPUT_FILE,TrainingTrajectoryFun,FILE_MATERIAL_DATA,DATALOC)
%%
% =========================================================================
% FEvectorized_GEN3 — Generate FE training snapshots for parametric loads
% =========================================================================
% PURPOSE
%   Vectorized driver to run the FE solver over a set of parametric load
%   states (the “training trajectory”) and store common data for downstream
%   OFFLINE/HROM stages.
%
% SIGNATURE
%   FEvectorized_GEN3(NAME_INPUT_FILE, TrainingTrajectoryFun, FILE_MATERIAL_DATA, DATALOC)
%
% INPUTS
%   NAME_INPUT_FILE        (char)  Name of FE input pack (e.g., 'INPUTS_FE_LOC').
%   TrainingTrajectoryFun  (char)  Function name that returns DATAcommon struct
%                                   with at least:
%                                     • INPUTS_PARAMETERS : array of load states
%                                     • NameParamStudy    : base label for runs
%                                     • SMALL_STRAIN_KINEMATICS (optional, default=1)
%   FILE_MATERIAL_DATA     (char)  Material file (e.g., 'PlasticMaterialDATA').
%   DATALOC                (struct) Misc runtime/config fields; augmented inside.
%
% DEFAULTS (when arguments are missing)
%   nargin==0 → NAME_INPUT_FILE='INPUTS_1elem'; TrainingTrajectoryFun='TrainingTrajectory';
%               FILE_MATERIAL_DATA='ElasticMaterialDATA' (or 'ElasticMaterialDATA1mat').
%   nargin==1 → TrainingTrajectoryFun='TrainingTrajectory';
%               FILE_MATERIAL_DATA='ElasticMaterialDATA'.
%   nargin==2 → FILE_MATERIAL_DATA='ElasticMaterialDATA'.
%
% WHAT THE FUNCTION DOES
%   1) Ensures code paths (env var MATLAB_CODES_FEHROM) if core FE functions
%      like 'Quadratic3N' are not on the path.
%   2) Logs console output to TRAINING.txt (diary) for reproducibility.
%   3) Builds DATAcommon by calling TrainingTrajectoryFun(DATALOC) and
%      enforces defaults via DefaultField(...).
%   4) Copies relevant flags to DATALOC (e.g., SMALL_STRAIN_KINEMATICS) and
%      sets DATALOC.StoreInitialStiffnessMatrix=1.
%   5) Loops over DATAcommon.INPUTS_PARAMETERS:
%        • Sets DATALOC.PARAM and DATALOC.ILOADSTATES.
%        • Tags each run with DATALOC.NameParamStudy = [NameParamStudy '_param_i'].
%        • Calls FECODElargeDYN(NAME_INPUT_FILE, DATALOC) to run the FE solve.
%        • On the first load state, stores common FE objects to
%          DATA.FE_VARIABLES_NAMEstore: MESH, MATPRO, OPERFE, DISP_CONDITIONS, OTHER_output.
%   6) Closes the diary.
%
% OUTPUTS / SIDE EFFECTS
%   • Files: TRAINING.txt (console log).
%   • MAT file (name in DATA.FE_VARIABLES_NAMEstore) containing common FE data:
%       MESH, MATPRO, OPERFE, DISP_CONDITIONS, OTHER_output.
%   • No function return. Results are consumed by offline builders / postprocs.
%
% KEY DEPENDENCIES
%   - FECODElargeDYN        : main FE solver.
%   - DefaultField          : utility for defaulting struct fields.
%   - TrainingTrajectoryFun : user-supplied function returning DATAcommon.
%   - Environment variable  : MATLAB_CODES_FEHROM must point to codebase root.
%
% PRACTICAL NOTES
%   • Set DATALOC.FILE_MATERIAL_DATA externally or via input; it is copied into
%     DATALOC internally for the FE solver.
%   • If RUN_PROBLEM is set to 0 (dev use), LocalPostProcess(NameParamStudy)
%     is called instead of running the FE solver.
%   • SMALL_STRAIN_KINEMATICS is enforced to default=1 unless supplied by
%     DATAcommon.
%
% USAGE EXAMPLE
%   FEvectorized_GEN3('INPUTS_FE_LOC','BCfe_plateT','PlasticMaterialDATA',DATALOC);
%
% .mlx references: (none)
%
% Automatically commented by ChatGPT
% =========================================================================

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
for ILOADSTATES =1:length(DATAcommon.INPUTS_PARAMETERS)
    DATALOC.PARAM  = DATAcommon.INPUTS_PARAMETERS(ILOADSTATES) ;
    DATALOC.ILOADSTATES = ILOADSTATES; 
    DATALOC.NameParamStudy =  [DATAcommon.NameParamStudy,'_param_',num2str(ILOADSTATES)]  ;
    if RUN_PROBLEM == 1
        [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] = ...
            FECODElargeDYN(NAME_INPUT_FILE,DATALOC) ;
        if ILOADSTATES == 1
            save(DATA.FE_VARIABLES_NAMEstore,'MESH','MATPRO','OPERFE','DISP_CONDITIONS','OTHER_output') ;
        end
    else
        % SHOW SVD ERRORS OF EACH SNAPSHOT
        LocalPostProcess(DATALOC.NameParamStudy) ;
    end
end
diary off
