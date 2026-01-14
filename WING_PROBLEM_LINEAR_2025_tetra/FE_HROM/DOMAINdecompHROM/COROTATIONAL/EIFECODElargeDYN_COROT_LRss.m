function [MESH,MATPRO,OPERFE,DATA,TIME_COMPUT,DISP_CONDITIONS,OTHER_output] =  EIFECODElargeDYN_COROT_LRss(DATAFILE,PARAMETERS)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: EIFECODElargeDYN_COROT_LRss
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Main driver for solving dynamic or static nonlinear mechanical problems using the **Empirical
%   Interscale Finite Element Method (EIFEM)** in the **corotational setting** for **large rotations
%   and small (or large) strains**. It orchestrates all steps from data loading, initialization,
%   time stepping, and postprocessing with GiD.
%
%   This is an extension of `EIFECODElargeDYN`, adapted for the COROT_LRss regime.
%
% USAGE:
%   [MESH, MATPRO, OPERFE, DATA, TIME_COMPUT, DISP_CONDITIONS, OTHER_output] = ...
%        EIFECODElargeDYN_COROT_LRss(DATAFILE, PARAMETERS)
%
% INPUTS:
%   - DATAFILE     : String with the name of the user-defined input file (as a function handle or string).
%                    This file must return mesh, material properties, FE operators, boundary conditions, etc.
%   - PARAMETERS   : Structure with simulation options, flags (e.g. ONLY_PRINT_GID), and paths.
%
% OUTPUTS:
%   - MESH             : Structure with mesh connectivity, geometry, element types, etc.
%   - MATPRO           : Structure with material properties used in constitutive modeling.
%   - OPERFE           : Finite element operators (boolean matrices, basis function info, etc.).
%   - DATA             : Main control structure with solver settings, results, and mesh info.
%   - TIME_COMPUT      : Structure measuring time taken by key phases (input, time loop, etc.).
%   - DISP_CONDITIONS  : Structure with displacement BCs, projection operators (A, G), and prescribed motion.
%   - OTHER_output     : Miscellaneous outputs (force histories, postprocess variables, stiffness, etc.).
%
% FUNCTIONALITY:
%   - Reads or generates all inputs via a user-defined file (`DATAFILE`).
%   - Initializes all simulation variables, including rotation matrices `Qrot`.
%   - Depending on `DATA.ISDYNAMIC`, it runs either:
%       * `NONLINEARstaticLARGE_COROT_LRss` – Static simulation in corotational frame.
%       * `NONLINEARdynamicLARGE` – Dynamic simulation (not implemented here).
%   - Stores time and iteration-level performance in `TIME_COMPUT`.
%   - Automatically triggers GiD postprocessing using:
%       * `GidPostProcess_VECT_EIFEbubCOROT` or `GidPostProcess_VECT_EIFE`
%       * Optionally writes binary GiD results if `SAVE_AS_BINARY_FILE = 1`
%   - Handles additional postprocessing of unconverged steps and rigid-body motions.
%
% FEATURES:
%   - Full corotational support for large rotation/small strain problems.
%   - Compatible with SVD-based space-time decomposition for body and traction forces.
%   - Seamless integration with EIFEM-specific element-level operations and projection-based DOFs.
%   - Supports precomputed stiffness matrices, internal force projection, and reduced models.
%
% REFERENCES:
%   - Theory and implementation details:
%     /PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
%   - Code structure inspired by:
%     * EIFECODElargeDYN.m
%     * EIFECODElargeDYN_COROT2.m
%     * Examples in /TESTING_PROBLEMS_FEHROM/...
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, UPC – CIMNE  
%   Date: 05-Feb-2024, Balmes 185, Barcelona
%   Comments created by ChatGPT on 12-May-2025
%
% DEPENDENCIES:
%   - DATAFILE (user-defined input generator)
%   - INITIALIZATIONvar
%   - ComputeKiniCheck, ComputeCini_B_COROT
%   - NONLINEARstaticLARGE_COROT_LRss
%   - GidPostProcess_VECT_EIFEbubCOROT / GidPostProcess_VECT_EIFE
%   - GidPostProcessLARGE_iterations (for unconverged steps)
%   - DefaultField (utility)
%   - GIDres2bin (optional binary output)
%
% ---------------------------------------------------------------------------------------------------




% EMPIRICAL INTERSCALE FINITE ELEMENT CODE,  dynamic regime, co-rotational,
% large strains 
% Copy of EIFECODElargeDYN, and then of  EIFECODElargeDYN_COROT.m, and finally
% of EIFECODElargeDYN_COROT2
% JAHO, 05-02-2024, Balmes 185,  Barcelona
% Theory developed in
% /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
% .m 
if nargin ==0
    load('tmp3.mat')
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

if PARAMETERS.ONLY_PRINT_GID == 0
    
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
    PARAMETERS = DefaultField(PARAMETERS,'OnlyComputeMassAndStiffnessMatricesSmallStrains',0) ; % 19-Feb-2022.
    OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[]) ; % 25-Feb-2022.
    PARAMETERS = DefaultField(PARAMETERS,'INTERNAL_FORCES_USING_precomputed_Kstiff',0) ; % 24-May-2022. Avoid computing stresses in linear problems
    DATA = DefaultField(DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff',PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff) ; % 24-May-2022. Avoid computing stresses in linear problems
    DATA = DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
    
    % Some variants require computing the initial stiffness matrix, among
    % other variables
     [OTHER_output,OPERFE]= ComputeKiniCheck(DATA,PARAMETERS,OTHER_output,OPERFE,MATPRO) ; 
     
    
    disp('*************************************')
    disp('INITIALIZATIONS')
    disp('*************************************')
    % Variables at time n (and previous time steps if required. Initialization of
    % snapshots )
    [DATA,VAR,SNAP] = INITIALIZATIONvar(DATA,INICOND)   ;
    
    [ OPERFE]= ComputeCini_B_COROT(DATA,OPERFE,MATPRO,VAR) ; 

    %
    
    %d = VAR_n.DISP(:,istep-1) ;
    TIME_COMPUT.TIME_STEP_LOOP = tic ;
    PARAMETERS = DefaultField(PARAMETERS,'ONLY_PRINT_GID',0) ;
    if PARAMETERS.ONLY_PRINT_GID == 0
        if DATA.ISDYNAMIC == 0
            [DATA,CONVERGED,QrotTIME]=NONLINEARstaticLARGE_COROT_LRss(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
        else
            error('Option not implemented yet')
            [DATA,CONVERGED]=NONLINEARdynamicLARGE(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO) ;
        end
    else
        disp('...Using results already calculated..')
    end
    TIME_COMPUT.TIME_STEP_LOOP = toc(TIME_COMPUT.TIME_STEP_LOOP) ;
    OTHER_output.TIME_COMPUT = TIME_COMPUT;
    OTHER_output.DATA = DATA ; % 25-May-2022 ************************
    disp(['...done in (',num2str(TIME_COMPUT.TIME_STEP_LOOP ),' s)'])
    
    
else
    disp('RETRIEVING INFORMATION FOR PRINTING GID FILES')
    FOLDER = [cd,filesep,'SNAPSHOTS',filesep]  ;
    FE_VARIABLES_NAMEstore = [FOLDER,filesep,PARAMETERS.NameParamStudy,'_FEoper','.mat'] ;
    
    load(FE_VARIABLES_NAMEstore,'MESH','OTHER_output') ;
    DATA = OTHER_output.DATA;
    DISP_CONDITIONS = [] ;
    MATPRO = [] ; OPERFE = [] ; TIME_COMPUT = [] ; DISP_CONDITIONS = [] ;
end





disp('Postprocess with GID')
% POST-PROCESS WITH GID

DATAinpGID = [] ;
PARAMETERS = DefaultField(PARAMETERS,'PRINT_RIGID_BODY',1 ) ;
if PARAMETERS.PRINT_RIGID_BODY == 1 && ~isempty(DISP_CONDITIONS)
    DATAinpGID.ADDITION_VARIABLE.DISP = DISP_CONDITIONS.RIGID_BODY_MOTION ;
end
 for icluster = 1:length(DATA.STORE.NSTEPS_CLUSTER)
    if  isempty(DATA.FOLLOWER_LOADS) ||  isempty(DATA.FOLLOWER_LOADS.HYDROSTATIC.NAME_MESH_WATERLINE)
        %   GidPostProcessLARGE(MESH,DATA,icluster,DATAinpGID);
        disp('Printing GID files ...')
        DATA =DefaultField(DATA,'MESHextended',[]);
        if ~isempty(DATA.MESHextended)
       OTHER_output.strainCOARSE_history =   ...
           GidPostProcess_VECT_EIFEbubCOROT(MESH,DATA,icluster,DATAinpGID,OTHER_output.PROPMAT,PARAMETERS,QrotTIME,QrotTIME{1},...
           OPERFE.LboolCall);
        else
            error('Option  not implemented')
         OTHER_output.strainCOARSE_history =    GidPostProcess_VECT_EIFE(MESH,DATA,icluster,DATAinpGID,OTHER_output.PROPMAT,PARAMETERS);
         % This output is temporary....30-Oct-2023
        end
        return
        
    else
        error('Option not implemented')
        %    GidPostProcessLARGE_waterline(MESH,DATA,icluster,DATAinpGID);
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


