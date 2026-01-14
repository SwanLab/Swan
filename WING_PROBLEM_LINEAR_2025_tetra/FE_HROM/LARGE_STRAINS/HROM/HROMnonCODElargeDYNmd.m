function  [BasisUall,RECONS_PK2stress,OTHER_output ]=  HROMnonCODElargeDYNmd(DATAFILE,PARAMETERS)
%--------------------------------------------------------------------------
% HROMnonCODElargeDYNmd
% ---------------------
% DESCRIPTION (updated documentation)
%   Executes a nonlinear Hyper-Reduced Order Model (HROM) simulation for
%   geometrically nonlinear structural problems with potentially multiple
%   generalized coordinates. This version extends
%   HROMnonCODElargeDYNfst by supporting multi-coordinate manifolds,
%   adaptive hyperreduction strategies, and optional GID post-processing.
%
% PURPOSE
%   1) Read and assemble all reduced-order input data (basis, mesh, ECM data,
%      material properties, and solver parameters).
%   2) Initialize variables and precompute tangent operators (C_ini) when
%      required.
%   3) Solve the nonlinear equilibrium problem (static or dynamic) projected
%      onto the tangent manifold using the reduced basis.
%   4) Reconstruct PK2 stresses and, optionally, internal variables for
%      postprocessing and visualization.
%
% KEY FEATURES
%   • Multi-coordinate tangent-manifold formulation.
%   • Hyperreduction via element-free ECM data structures.
%   • Optional adaptive weights (SAW-ECM strategy).
%   • Modular initialization: reads data from an input generator (DATAFILE).
%   • Optional GiD export with reduced mesh visualization and stress fields.
%
% INPUTS
%   DATAFILE   : Function handle returning model data:
%                [BasisUall, BasisStwo, ECMdata, MESH, MATPRO, OPERHROM,
%                 Fbody, Ftrac, DISP_CONDITIONS, INICOND, DATAHROM, OTHER_output]
%   PARAMETERS : Structure containing solver, integration, and output options.
%
% OUTPUTS
%   BasisUall        : Reduced displacement basis used for reconstruction.
%   RECONS_PK2stress : Reconstructed 2nd Piola–Kirchhoff stress field.
%   OTHER_output     : Structure aggregating auxiliary data:
%                      - TIME_COMPUT: timing information
%                      - DISP_CONDITIONS, OPERHROM, DATAHROM
%                      - optional ECM cluster weights, mesh info, etc.
%
% EXECUTION FLOW
%   1) INPUT READING
%      - Calls DATAFILE(PARAMETERS) to obtain reduced model data.
%      - Initializes OPERHROM and DATAHROM optional fields using DefaultField.
%
%   2) INITIALIZATION
%      - Builds initial variables VAR and SNAP via INITIALIZATIONvar.
%      - Computes initial stiffness (C_ini) for linear-elastic portions if
%        DATAHROM.CECM_ONLY_FOR_NONLINEAR_STRESSES = 1.
%
%   3) SOLVER STAGE
%      - If ISDYNAMIC = 0 → calls NONLINEARstaticLARGE_manifoldMD.
%      - If ISDYNAMIC = 1 → (dynamic solver placeholder).
%      - Records computation time.
%
%   4) RECONSTRUCTION AND POSTPROCESSING
%      - Builds DATAinpGID structure for GiD export.
%      - Optionally reconstructs PK2 stresses via PK2stressReconstruction.
%      - Optionally prints GiD results (GidPostProcessLARGE or
%        GidPostProcessLARGE_waterline), and can convert .res → .bin.
%
% NOTES / IMPLEMENTATION DETAILS
%   - Compatible with adaptive ECM weighting schemes (OPERHROM.wSTs_cluster).
%   - If Show_HROM_MESH=1, adds ECM visualization to GiD output.
%   - Optional binary export controlled by DATAHROM.PRINT.SAVE_AS_BINARY_FILE.
%   - Time-tracking fields in OTHER_output.TIME_COMPUT:
%       INPUTS, TIME_STEP_LOOP (seconds).
%
% DEPENDENCIES
%   - INITIALIZATIONvar, ComputeCini_B, NONLINEARstaticLARGE_manifoldMD,
%     PK2stressReconstruction, GidPostProcessLARGE(_waterline),
%     GIDres2bin, DefaultField, etc.
%
% REFERENCES
%   - Hernández Ortega, J.A. (2025). "Machine Learning Techniques in
%     Structural Analysis" (Sections 9–11): tangent-space projection,
%     nonlinear parametrization τ(q), and hyperreduction.
%   - Eq. (176–236): projection and stress-reconstruction procedures.
%   - See implementation examples in:
%       /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/
%
% VERSION HISTORY
%   - Derived from HROMnonCODElargeDYNfst
%   - Date modification: 13th August 2025, Wednesday. Molinos Marfagones, Cartagena.
%   - Author: Joaquín A. Hernández Ortega (JAHO)
%     Original creation: 2 July 2025, Honest Greens Pedralbes, Barcelona.
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------


if nargin ==0
    load('tmp.mat')
end
TIME_COMPUT=[] ;
format long g

% READING AND CONSTRUCTING INPUT DATA FOR THE FE CODE
% -----------------------------------------------------------
% Prototypical file, see
%/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/HROM/BEAM_2inpvar_SMALL_DISP/INPUTS_HROM_2param.m
disp('*************************************')
disp('Reading and constructing HROM input data')
disp('*************************************')
TIME_COMPUT.INPUTS = tic ;
[BasisUall,BasisStwo,ECMdata,...
    MESH,MATPRO,OPERHROM,Fbody,Ftrac,DISP_CONDITIONS,INICOND,DATAHROM,OTHER_output]  = ...
    feval(DATAFILE,PARAMETERS) ;

OPERHROM = DefaultField(OPERHROM,'HYDRO',[]) ;
OPERHROM = DefaultField(OPERHROM,'DOFm',[]) ;  % 27-Jan-2024


DATAHROM = DefaultField(DATAHROM,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ;

OTHER_output.DISP_CONDITIONS =DISP_CONDITIONS ;
OTHER_output.OPERHROM = OPERHROM;

OTHER_output.DATAHROM = DATAHROM ;
%OTHER_output.Ftrac = Ftrac; 

TIME_COMPUT.INPUTS = toc(TIME_COMPUT.INPUTS) ;
disp(['...done in (',num2str(TIME_COMPUT.INPUTS ),' s)'])
disp('*************************************')

% ------------------------------------------------------------------------------------------------------------

disp('*************************************')
disp('INITIALIZATIONS')
disp('*************************************')
% Variables at time n (and previous time steps if required. Initialization of
% snapshots )
[DATAHROM,VAR,SNAP] = INITIALIZATIONvar(DATAHROM,INICOND)   ;
% 24th-Aug-2025, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/08_PLAST_adapECM.mlx
% Computation initial elasticity matrix for cases in which
% DATAHROM.CECM_ONLY_FOR_NONLINEAR_STRESSES = 1
VAR_TMP=VAR; 
  VAR_TMP.DISP = zeros(size(OPERHROM.Bst,2),1) ; 
[ OPERHROM]= ComputeCini_B(DATAHROM,OPERHROM,MATPRO,VAR_TMP) ; 
 %
PARAMETERS = DefaultField(PARAMETERS,'INTERNAL_FORCES_USING_precomputed_Kstiff',0) ; % 24-May-2022. Avoid computing stresses in linear problems
DATAHROM = DefaultField(DATAHROM,'INTERNAL_FORCES_USING_precomputed_Kstiff',PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff) ; % 24-May-2022. Avoid computing stresses in linear problems


% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/09_SAW_ECMplast.mlx
% This is for adaptive WEIGHTS using SAW-ECM
if isstruct(OPERHROM.wSTs)
    OPERHROM.wSTs_cluster = OPERHROM.wSTs;  
else
    OPERHROM.wSTs_cluster  = [] ; 
end

%d = VAR_n.DISP(:,istep-1) ;
TIME_COMPUT.TIME_STEP_LOOP = tic ;

PARAMETERS = DefaultField(PARAMETERS,'ONLY_PRINT_GID',0) ;
if PARAMETERS.ONLY_PRINT_GID == 0
    if DATAHROM.ISDYNAMIC == 0
        [DATAHROM,CONVERGED]=NONLINEARstaticLARGE_manifoldMD(DATAHROM,DISP_CONDITIONS,VAR,OPERHROM,SNAP,Fbody,Ftrac,MATPRO) ;
    else
        error('Option not implemented yet')
        [DATAHROM,CONVERGED]=NONLINEARdynamicLARGE(DATAHROM,DISP_CONDITIONS,VAR,OPERHROM,SNAP,Fbody,Ftrac,MATPRO) ;
    end
else
    disp('...Using results already calculated..')
end
TIME_COMPUT.TIME_STEP_LOOP = toc(TIME_COMPUT.TIME_STEP_LOOP) ;
OTHER_output.TIME_COMPUT = TIME_COMPUT;
disp(['...done in (',num2str(TIME_COMPUT.TIME_STEP_LOOP ),' s)'])


% POST-PROCESS WITH GID
% RECONSTRUCTION DISPLACEMENTS
% -----------------------------
DATAinpGID = [];
DATAinpGID.OPERreconstr.DISP.BASIS =  BasisUall ;
DATAinpGID.OPERreconstr.DISP.coeff =  1 ;

% Reconstruction of PK2 stresses
 % % Reconstruction of internal variables (if any)  %pending----
% if ~isempty(DATAHROM.ListFieldInternalVariables)
%     [OTHER_output,DATAinpGID,RECONS_PK2stress] = InternalVarReconstruction(DATAinpGID,OTHER_output,ECMdata,DATAHROM,BasisStwo) ;
% end
if ~isempty(BasisStwo)
[OTHER_output,DATAinpGID,RECONS_PK2stress] = PK2stressReconstruction(DATAinpGID,OTHER_output,ECMdata,DATAHROM,BasisStwo) ;
else
    RECONS_PK2stress = [] ; 
end


PARAMETERS = DefaultField(PARAMETERS,'PRINT_RIGID_BODY',1 ) ;
if PARAMETERS.PRINT_RIGID_BODY == 1
    DATAinpGID.ADDITION_VARIABLE.DISP = DISP_CONDITIONS.RIGID_BODY_MOTION ;
end
PARAMETERS = DefaultField(PARAMETERS,'Show_HROM_MESH',0 ) ;  % Show HROM mesh
if PARAMETERS.Show_HROM_MESH ==1
    DATAinpGID.ECMdata = ECMdata ;
    if isfield(OTHER_output,'ECMdata_press')
        DATAinpGID.ECMdata_press = OTHER_output.ECMdata_press ;
    end
    DATAinpGID.IS_HROM_SIMULATION = 1;
else
    DATAinpGID.IS_HROM_SIMULATION  =0;
end

PARAMETERS = DefaultField(PARAMETERS,'PRINT_POSTPROCESS_GID',0) ;

if PARAMETERS.PRINT_POSTPROCESS_GID == 1
    
    for icluster = 1:length(DATAHROM.STORE.NSTEPS_CLUSTER)
        % GidPostProcessLARGE(MESH,DATAHROM,icluster,DATAinpGID);
        
        if  isempty(DATAHROM.FOLLOWER_LOADS) ||  isempty(DATAHROM.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE)
            GidPostProcessLARGE(MESH,DATAHROM,icluster,DATAinpGID);
        else
            GidPostProcessLARGE_waterline(MESH,DATAHROM,icluster,DATAinpGID);
        end
        
        
    end
    
    
    
    DATAHROM.SNAP_ITER = [] ;
    if iscell(DATAHROM.PRINT.NAME_FILE_MSH)
        fff= DATAHROM.PRINT.NAME_FILE_MSH{1} ;
    else
        fff= DATAHROM.PRINT.NAME_FILE_MSH ;
    end
    BASE_FOLDER = fileparts(fff) ;
    DATAHROM.PRINT = DefaultField(DATAHROM.PRINT,'BASE_FOLDER',BASE_FOLDER)  ;
    DATAHROM.PRINT = DefaultField(DATAHROM.PRINT,'SAVE_AS_BINARY_FILE',0)  ;
    if DATAHROM.PRINT.SAVE_AS_BINARY_FILE == 1
        GIDres2bin(DATAHROM.PRINT.BASE_FOLDER) ;
    end
    
end


