function [CECMoutput,MSG,Wdom] = CECM_multi(DATAIN,DATAROM,BASES,DATA_REFMESH)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/MultiscaleHROM/READMEmulti.mlx
% as well as the reduced Bmatrix



if nargin == 0
    EXECUTABLE_FOLDER = getenv('MATLAB_CODES_FEHROM') ;
    
    if exist('ContinuousECMgen') ==0
        addpath(genpath(EXECUTABLE_FOLDER)) ;
    end
    
    cd('/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/MultiscaleHROM/') ;
    %load('CheckCECM2.mat')  % Periodic, quadratic
    %  load('CheckCECM_b.mat')  % Non-periodic, quadratic
    %     load('CheckCECM_c.mat')  % Non-periodic, linear 3x3x3 Gauss rule
    %     load('CheckCECM_d.mat')  % periodic, quadratic 5x5x5 Gauss rule
        load('CheckCECM_g.mat')  % I-shaped, quadratic, periodic, 27-points rule
   % load('CheckCECM_h.mat')   % I-shaped, quadratic, periodic, 64-points rule
   %   load('CheckCECM_i.mat')   % I-shaped, quadratic, periodic, 125-points rule, coarse
   %    load('CheckCECM_j.mat')   % I-shaped, linear, periodic, 8-points rule 
    DATAIN.CUBATURE.DATA_ECM.TOL_SVD_A = 1e-2;
    DATAIN.CUBATURE.DATA_ECM.EXCLUDE_ELEMENTS_TRIGGERED_NONCONVERGENCE_ITERATION_IN_CECM = 0;
    DATAIN.CUBATURE.DATA_ECM.SECOND_STAGE_ITERATIONS = 20;
end
% Continuous ECM
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'DATA_ECM',[]) ;
DATA_ECM = DATAIN.CUBATURE.DATA_ECM ;

% Default values
DATA_ECM = DefaultField(DATA_ECM,'TOL_SVD_A',DATAIN.CUBATURE.TOL_LOC_InternalForces) ;
DATA_ECM = DefaultField(DATA_ECM,'UsePartitionedRandomizedAlgorithm',0) ;
DATA_ECM = DefaultField(DATA_ECM,'ListElementsExclude_fromGID',[]) ;
DATA_ECM = DefaultField(DATA_ECM,'TOL_NewtonRaphson_EliminationPoints',1e-8) ;
DATA_ECM = DefaultField(DATA_ECM,'MaxIterationsNR_ElimPoints',40) ;
DATA_ECM = DefaultField(DATA_ECM,'SECOND_STAGE_ITERATIONS',20) ;
DATA_ECM = DefaultField(DATA_ECM,'Include2ndStageIterations_PlotEvolutionWeights',1) ;
DATA_ECM = DefaultField(DATA_ECM,'PLOT_INTERNAL_FORCE_SNAPSHOTS',1) ;
DATA_ECM = DefaultField(DATA_ECM,'PLOT_INTERNAL_MODES',1) ;

DATA_ECM = DefaultField(DATA_ECM,'SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY',0) ;

DATA_from_MAIN.SHOW_ALSO_DECM_ITERATIVE_PROCESS_GRAPHICALLY = 0 ;


disp('****************************+')
disp('Continuous Empirical Cubature Method ')
disp('*************************************')
load(DATAIN.NAME_WS_MODES,'CgloDOM','Wdom','Bdom')

MSG = {} ;
BdomRED = Bdom*DATAROM.BasisUdef ;

% ------------------------------
% Matrix of internal forces (integrand functions)
% -----------------------------
[A,MSG] = GetInternalForceModes(BASES,DATAIN,MSG,BdomRED,Wdom,DATAROM,Bdom) ;
% Plotting the matrix of integrand snapshots
if DATA_ECM.PLOT_INTERNAL_FORCE_SNAPSHOTS == 1
    NameFileMesh = [DATAIN.NAME_WS_MODES] ; % DATAIN.NAME_project ;
    DATAIN.LABEL_NAME_FILE = '_snapFINT' ;
    DATAIN =   GidPostProcessModesFINT_LOCg(DATA_REFMESH.COOR,DATA_REFMESH.CN,...
        DATA_REFMESH.TypeElement,A,DATA_REFMESH.posgp,NameFileMesh,DATAIN);
end
% Exact integral
DATA.ExactIntegral = A'*Wdom ;

% FINITE ELEMENT MESH  ---
% --------------------%
xNODES= DATA_REFMESH.COOR' ;
ndim = size(DATA_REFMESH.COOR,2) ;
xNODES = xNODES(:) ;
xFE= DATA_REFMESH.Nst*xNODES;
xFE =reshape(xFE,ndim,[])' ;

% GID elements excluded from the set of candidates (empty by default)
DATA_ECM = DefaultField(DATA_ECM,'ListElementsExclude_fromGID',[]) ;
if ~isempty(DATA_ECM.ListElementsExclude_fromGID)
    ListElementsToExclude = load(DATA_ECM.ListElementsExclude_fromGID) ;
    DATA_ECM.ListElementsToExclude = ListElementsToExclude(:,1) ;
else
    DATA_ECM.ListElementsToExclude = [] ;
end
% INVOKING THE CECM
% --------------------
DATA_REFMESH.ngausE = size(DATA_REFMESH.posgp,2) ;
DATAloc.ExactIntegral = A'*Wdom ;
DATAloc.BASE_FOLDER  = [DATAIN.LOCAL_FOLDER,filesep,'MODES',filesep] ;
[CECMoutput,DATA_AUX]= ContinuousECMgen(A,xFE,Wdom,DATAloc,DATA_REFMESH,DATA_ECM) ;

% PLOTTING THE MODES (IN GID)
if DATA_ECM.PLOT_INTERNAL_MODES == 1
    NameFileMesh = [DATAIN.NAME_WS_MODES] ; % DATAIN.NAME_project ;
    DATAIN.LABEL_NAME_FILE = '_MODESfint' ;
    DATAIN =   GidPostProcessModesFINT_LOCg(DATA_REFMESH.COOR,DATA_REFMESH.CN,...
        DATA_REFMESH.TypeElement,DATA_AUX.VAR_SMOOTH_FE.BasisIntegrand,DATA_REFMESH.posgp,NameFileMesh,DATAIN);
end

disp('----------------------------')
disp('Stiffness matrix DECM ')
disp('---------------------------')
nstrain = size(BdomRED,1)/length(Wdom);
DATAlocS = [] ;
[Kstiff_DECM,~,~] = StiffReduced(DATA_AUX.indexPoints_DECM, nstrain,CgloDOM,DATAIN,BdomRED,CECMoutput.wDECM,Wdom,DATAlocS) ;
% Comparison with the one computed used the whole set of Gauss points
Kstiff_all = DATAROM.KdomRED ;
DIFFERENCE = norm((Kstiff_DECM-Kstiff_all),'fro')/norm(Kstiff_all,'fro')*100 ;
MSGloc1 = ['Difference between HROM-DECM and FE Stiff Matrix = ',num2str(DIFFERENCE),' %'] ;
disp(MSGloc1)
MSG{end+1} = [MSGloc1,'(per cent)'] ;
disp('---------------------------')
disp('----------------------------')
disp('Stiffness matrix CECM ')
disp('---------------------------')
Kstiff_CECM = StiffReducedCECM(CECMoutput, nstrain,CgloDOM,DATAIN,BdomRED,Wdom,DATA_AUX) ;
DIFFERENCE = norm((Kstiff_CECM-Kstiff_all),'fro')/norm(Kstiff_all,'fro')*100 ;
MSGloc1 = ['Difference between HROM-CECM and FE Stiff Matrix = ',num2str(DIFFERENCE),' %'] ;
disp(MSGloc1)
MSG{end+1} = [MSGloc1,'(per cent)'] ;


CECMoutput.nstrain = nstrain ;

%