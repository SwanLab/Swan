function [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] =  FECODElargeDYN(DATAFILE,PARAMETERS)
%%
% =========================================================================
% FECODElargeDYN — Finite Element solver (large strains, static/dynamic)
% =========================================================================
% PURPOSE
%   General FE driver that assembles and solves nonlinear problems at large
%   strains in static or dynamic regime, with optional precomputations for
%   small-strain stiffness and hooks for HROM/ECM workflows.
%
% SIGNATURE
%   [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,CONVERGED,DISP_CONDITIONS,OTHER_output] = ...
%       FECODElargeDYN(DATAFILE, PARAMETERS)
%
% INPUTS
%   DATAFILE    (char)   Name of a function returning the FE input pack:
%                        [MESH,MATPRO,OPERFE,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATA,OTHER_output].
%   PARAMETERS  (struct) Runtime controls; missing fields are defaulted via DefaultField.
%                        Recognized fields include (not exhaustive):
%                          • OnlyComputeMassAndStiffnessMatricesSmallStrains (logical)
%                          • INTERNAL_FORCES_USING_precomputed_Kstiff        (logical)
%                          • StoreInitialStiffnessMatrix                      (logical)
%                          • ONLY_PRINT_GID                                   (logical)
%                          • PRINT_RIGID_BODY                                 (logical)
%
% OUTPUTS
%   MESH, MATPRO, OPERFE, DATA : Updated FE data structures after solve.
%   TIME_COMPUT                : Timers for input assembly and time-step loop.
%   CONVERGED                  : Convergence info/status returned by the solver.
%   DISP_CONDITIONS            : Dirichlet data (incl. rigid-body motion if requested).
%   OTHER_output               : Aux data (forces, initial matrices, timing, etc.).
%
% WHAT THE FUNCTION DOES
%   1) Read/construct inputs by calling feval(DATAFILE, PARAMETERS).
%   2) Default missing flags (e.g., SMALL_STRAIN_KINEMATICS, PRECOMPUTE_ELASTIC_STIFFNESS_MATRIX,
%      CECM_ONLY_FOR_NONLINEAR_STRESSES, INTERNAL_FORCES_USING_precomputed_Kstiff, etc.).
%   3) Initialize time-history variables and snapshots via INITIALIZATIONvar.
%   4) Precompute/check initial stiffness and related data via ComputeKiniCheck_VAR
%      (variant compatible with recent EIFECODElargeDYN logic).
%   5) Solve:
%        • If DATA.ISDYNAMIC == 0 → NONLINEARstaticLARGE(...)
%        • Else                    → NONLINEARdynamicLARGE(...)
%      (skips solve if PARAMETERS.ONLY_PRINT_GID == 1 and uses existing results).
%   6) Postprocess to GiD for each stored cluster:
%        • GidPostProcessLARGE or GidPostProcessLARGE_waterline
%        • Optionally print rigid-body motion (PARAMETERS.PRINT_RIGID_BODY)
%        • If non-converged iterations exist → GidPostProcessLARGE_iterations
%   7) Optional binary export of GiD results via GIDres2bin when requested.
%
% TIMERS
%   TIME_COMPUT.INPUTS         : Input assembly time.
%   TIME_COMPUT.TIME_STEP_LOOP : Nonlinear solve loop (static/dynamic).
%
% KEY DEPENDENCIES
%   DATAFILE (problem pack), DefaultField, INITIALIZATIONvar, ComputeKiniCheck_VAR,
%   NONLINEARstaticLARGE, NONLINEARdynamicLARGE, GidPostProcessLARGE(_waterline/_iterations),
%   GIDres2bin.
%
% PRACTICAL NOTES
%   • PARAMETERS and DATA flags are harmonized so linear cases can reuse
%     precomputed K (INTERNAL_FORCES_USING_precomputed_Kstiff) to avoid stress loops.
%   • For reproducibility, upstream drivers typically enable diary logging
%     (e.g., TRAINING.txt) outside this function.
%   • BASE_FOLDER for GiD output is inferred from DATA.PRINT.NAME_FILE_MSH.
%
% .mlx references: (none)
%
% Dates/places (original header retained for traceability)
%   Original: JAHO, 26-Nov-2020
%
% Written by Joaquín A. Hernández, UPC/CIMNE
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

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
DATA = DefaultField(DATA,'PRECOMPUTE_ELASTIC_STIFFNESS_MATRIX',0) ; % 9-Feb-2022.
DATA = DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; %

PARAMETERS = DefaultField(PARAMETERS,'OnlyComputeMassAndStiffnessMatricesSmallStrains',0) ; % 19-Feb-2022.
OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[]) ; % 25-Feb-2022.
OPERFE = DefaultField(OPERFE,'DOFm',[]) ; % 27-Feb-2024.

PARAMETERS = DefaultField(PARAMETERS,'INTERNAL_FORCES_USING_precomputed_Kstiff',0) ; % 24-May-2022. Avoid computing stresses in linear problems
DATA = DefaultField(DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff',PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff) ; % 24-May-2022. Avoid computing stresses in linear problems
PARAMETERS = DefaultField(PARAMETERS,'StoreInitialStiffnessMatrix',0) ; %  17-May-2024



%
%
% if (DATA.SMALL_STRAIN_KINEMATICS == 1 && DATA.PRECOMPUTE_ELASTIC_STIFFNESS_MATRIX  == 1) || ...
%         PARAMETERS.OnlyComputeMassAndStiffnessMatricesSmallStrains == 1 || PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff == 1 ...
%     || PARAMETERS.StoreInitialStiffnessMatrix  == 1
%     % Determine Stiffness Matrix, small strains
%   %  d  = zeros(DATA.MESH.ndof,1) ;
%     VAR.DISP =zeros(DATA.MESH.ndof,1) ;VARint_n =[] ;
%     [~,celastST,FgradST,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
%     OTHER_output.K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST,DATA) ;
%     if PARAMETERS.OnlyComputeMassAndStiffnessMatricesSmallStrains == 1
%         CONVERGED = 0 ;
%         return
%     end
%     if PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff == 1
%       OPERFE.KinternalFORCES_given = OTHER_output.K ;
%     end
% end
% 
disp('*************************************')
disp('INITIALIZATIONS')
disp('*************************************')
% Variables at time n (and previous time steps if required. Initialization of
% snapshots )
[DATA,VAR,SNAP] = INITIALIZATIONvar(DATA,INICOND)   ;

% 31-MAY-2O24, DONE IN ANALOGY TO
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/EIFECODElargeDYN.m
% If something goes wrong, uncomment the above portion of code
% Some variants require computing the initial stiffness matrix, among
% other variables
 %[OTHER_output,OPERFE]= ComputeKiniCheck(DATA,PARAMETERS,OTHER_output,OPERFE,MATPRO) ;
 % Modification introduced on 16-08-2025. If something goes wrong, you can
 % return to the previous implementation FECODElargeDYN_backup_16_08_2025.m
 % 
[OTHER_output,OPERFE]= ComputeKiniCheck_VAR(DATA,PARAMETERS,OTHER_output,OPERFE,MATPRO,VAR) ;
%



% 

%d = VAR_n.DISP(:,istep-1) ;
TIME_COMPUT.TIME_STEP_LOOP = tic ;
PARAMETERS = DefaultField(PARAMETERS,'ONLY_PRINT_GID',0) ;
if PARAMETERS.ONLY_PRINT_GID == 0
    if DATA.ISDYNAMIC == 0
        [DATA,CONVERGED]=NONLINEARstaticLARGE(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
    else
        [DATA,CONVERGED]=NONLINEARdynamicLARGE(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
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
PARAMETERS = DefaultField(PARAMETERS,'PRINT_RIGID_BODY',1 ) ;
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


