function Training_MULTI(INPUT_PROBLEMS,DATAoffline) 
%--------------------------------------------------------------------------
% FUNCTION: Training_MULTI
%
% PURPOSE:
%   This function performs the full-field finite element simulations
%   required to collect training data (snapshots) for the multiscale
%   Empirical Interscale Finite Element Method (EIFEM). It loops over
%   multiple load or parameter states defined in the `INPUT_PROBLEMS`
%   function and stores the resulting internal force components, body
%   forces, and boundary tractions.
%
%   This step constitutes the first phase (training stage) in a reduced
%   order modeling (ROM) pipeline, enabling later computation of basis
%   matrices and hyperreduction operators.
%
% INPUTS:
%   - INPUT_PROBLEMS : Function handle (string) that defines:
%         * Mesh file name
%         * Material parameters
%         * Dirichlet/Neumann boundary conditions
%         * Parameterized load trajectories (DATALOC.INPUTS_PARAMETERS)
%         * Flags for small strain or large strain settings
%
%   - DATAoffline    : (Optional) Structure allowing customization of
%                      mesh name, labels, and printing settings. The field:
%                         * DATAoffline.NAMEMESH
%                      overrides the default mesh file in DATALOC.
%
% OUTPUT:
%   - Snapshots of displacement, internal variables, stresses, etc. saved to disk.
%   - First workspace file (`FE_VARIABLES_NAMEstore`) contains:
%         * MESH: Mesh structure used for simulation
%         * MATPRO: Material property structure
%         * OPERFE: Finite element operator structure
%         * DISP_CONDITIONS: Dirichlet BC structure
%         * OTHER_output.Fbody: Cell array with body forces per param. state
%         * OTHER_output.Ftrac: Cell array with Neumann loads per param. state
%
% NOTES:
%   - Results are written in the folder specified in `DATALOC.FE_VARIABLES_NAMEstore`.
%   - Execution trace is saved in 'TRAINING.txt'.
%
% RELATED FILES:
%   - See implementation notes in:
%         /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/
%         TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Created: 2023-01-30 — Last Updated: 2023-02-25
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

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
DATALOC =    feval(INPUT_PROBLEMS,DATAoffline) ;    % Recovering the variables defined inside function "INPUT_PROBLEMS" 
DATAoffline = DefaultField(DATAoffline,'NAMEMESH', DATALOC.NameFileMeshDATA); 
DATALOC.NameFileMeshDATA = DATAoffline.NAMEMESH ; 

%DATAcommon =    feval(TrainingTrajectoryFun) ;   


DATALOC.FileUsedToRunParametricStudy = [cd,filesep,mfilename] ; 
% LOOP OVER PARAMETRIC SPACE
Fbody = cell(1,size(DATALOC.INPUTS_PARAMETERS,2) ) ; 
Ftrac = cell(1,size(DATALOC.INPUTS_PARAMETERS,2) ) ; 

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
end
save(NAME_FE_WS_1st,'Fbody','Ftrac','-append') ;
diary off
