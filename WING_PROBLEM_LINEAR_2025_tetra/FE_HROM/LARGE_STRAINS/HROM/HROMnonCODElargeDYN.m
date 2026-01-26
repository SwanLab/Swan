function  [BasisUall,RECONS_PK2stress,OTHER_output ]=  HROMnonCODElargeDYN(DATAFILE,PARAMETERS)
%--------------------------------------------------------------------------
% function [BasisUall, RECONS_PK2stress, OTHER_output] = HROMnonCODElargeDYN(DATAFILE, PARAMETERS)
%
% PURPOSE:
%   Executes a nonlinear hyperreduced-order model (HROM) simulation for a
%   dynamical or quasistatic structural problem with geometric nonlinearities.
%   The model follows the tangent-space projection formalism detailed in
%   Section 9.2 of MLEARNstruct_1.pdf, in which the reduced displacement
%   vector is obtained as a nonlinear function τ(q) of the generalized coordinates q.
%
%   The internal force vector is evaluated using an element-free formulation,
%   while the reduced model equations are solved via Galerkin projection onto
%   the tangent space to the solution manifold. Hyperreduction is employed
%   to accelerate the assembly of the reduced internal force vector using a
%   subset of Gauss points.
%
%   The function also reconstructs the PK2 stress field and optionally exports
%   simulation results to GID format.
%
% INPUTS:
%   - DATAFILE     : Function handle returning the reduced basis matrices, mesh,
%                    material data, force terms, boundary and initial conditions.
%   - PARAMETERS   : Structure specifying user options for printing, integration,
%                    solver tolerances, etc.
%
% OUTPUTS:
%   - BasisUall        : Matrix with global displacement basis vectors used for reconstruction.
%   - RECONS_PK2stress : Reconstructed PK2 stresses across the full mesh.
%   - OTHER_output     : Structure with auxiliary outputs including:
%                         • time performance
%                         • data structures (MESH, ECMdata, OPERHROM)
%                         • snapshots, boundary conditions, and solver results.
%
% REFERENCES:
%   - J.A. Hernández Ortega, "Machine Learning Techniques in Structural Analysis",
%     Sections 9 and 11 (pp. 71–97): Nonlinear ROM theory, tangent-space projection,
%     and hyperreduction strategies.
%
%   - Eq. (176)–(193): Nonlinear parametrization of reduced displacements via τ(q)
%     and residual projection onto tangent space.
%
%   - Eq. (232)–(236): Hyperreduced residual and stress recovery procedures.
%
%   - See also implementation examples in folder:
%     /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, 2 July 2025, Honest Greens Pedralbes, Barcelona.
%
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
%
PARAMETERS = DefaultField(PARAMETERS,'INTERNAL_FORCES_USING_precomputed_Kstiff',0) ; % 24-May-2022. Avoid computing stresses in linear problems
DATAHROM = DefaultField(DATAHROM,'INTERNAL_FORCES_USING_precomputed_Kstiff',PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff) ; % 24-May-2022. Avoid computing stresses in linear problems





%d = VAR_n.DISP(:,istep-1) ;
TIME_COMPUT.TIME_STEP_LOOP = tic ;

PARAMETERS = DefaultField(PARAMETERS,'ONLY_PRINT_GID',0) ;
if PARAMETERS.ONLY_PRINT_GID == 0
    if DATAHROM.ISDYNAMIC == 0
        [DATAHROM,CONVERGED]=NONLINEARstaticLARGE_manifold(DATAHROM,DISP_CONDITIONS,VAR,OPERHROM,SNAP,Fbody,Ftrac,MATPRO) ;
    else
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
[OTHER_output,DATAinpGID,RECONS_PK2stress] = PK2stressReconstruction(DATAinpGID,OTHER_output,ECMdata,DATAHROM,BasisStwo) ;
% % Reconstruction of internal variables (if any)  %pending----
% if ~isempty(DATAHROM.ListFieldInternalVariables)
%     [OTHER_output,DATAinpGID,RECONS_PK2stress] = InternalVarReconstruction(DATAinpGID,OTHER_output,ECMdata,DATAHROM,BasisStwo) ;
% end



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


