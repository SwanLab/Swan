function CECM= ContinuousECM_multi(OPERFE,PhiDEF,BasisPone,DATA,DATAcommon,MESH,DATAoffline,MATPRO,PhiRB) ;
% CONTINUOUS EMPIRICAL CUBATURE METHOD
% Invoked from, for instance,
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/01_INNER_DOM/OFFLINE_STEPSfun.m
% See further info in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% JAHO, 14-FEB-2023

if nargin == 0
    load('tmp3.mat')
    %  DATAoffline.LOAD_FROM_MEMORY_ECM_POINTS_INTERNAL_FORCES = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matrix of reduced internal forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('----------------------------------------')
disp('HYPERREDUCTION INTERNAL FORCES')
disp('------------------------------------------')

DATA_ECM.ACTIVE =1;  % Enable the continuous empirical cubature method
DATA_ECM.TOL_SVD_A =DATAoffline.errorFINT; % Tolerance SVD data matrix
DATA_ECM.UsePartitionedRandomizedAlgorithm = 1 ;  % Use advanced SVD
DATA_ECM.TOL_NewtonRaphson_EliminationPoints = 1e-8 ;  % Tolerance for the Newton-Raphson algorithm
DATA_ECM.MaxIterationsNR_ElimPoints = 40 ;  % Number of iteration for each NR step
SECOND_STAGE_ITERATIONS  = 20 ;   % Number of steps for each weight elimination in the second stage

DATAoffline = DefaultField(DATAoffline,'CECM_SECOND_STAGE_ITERATIONS',SECOND_STAGE_ITERATIONS) ;
DATA_ECM.SECOND_STAGE_ITERATIONS = DATAoffline.CECM_SECOND_STAGE_ITERATIONS ;


DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',[]) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'INTERNAL_FORCES',0) ;
DATAoffline = DefaultField(DATAoffline,'NumberOfCECMpoints',[]) ;
DATAoffline = DefaultField(DATAoffline,'INITIAL_SET_FOR_ECM_INDEX_MATERIAL',[]) ;
DATAoffline = DefaultField(DATAoffline,'IncludeSingularValueBasisStressesPone',0) ; % 
DATA = DefaultField(DATA,'BasisSTRESS_SINGULAR_VALUES',ones(size(BasisPone,2),1)) ; % 

if DATAoffline.IncludeSingularValueBasisStressesPone == 0
    DATA.BasisSTRESS_SINGULAR_VALUES = ones(size(BasisPone,2),1) ; 
end


if ~isempty(DATAoffline.INITIAL_SET_FOR_ECM_INDEX_MATERIAL)
    % INITIAL SET FOR ECM CANDITATE POINTS EQUAL TO THE POINTS
    % CORRESPONDING TO A GIVEN MATERIAL
    InitialSetElementsInternalForces = cell(length(DATAoffline.INITIAL_SET_FOR_ECM_INDEX_MATERIAL),1) ; %MESH.MaterialType ;
    for imat = 1:length(DATAoffline.INITIAL_SET_FOR_ECM_INDEX_MATERIAL)
        InitialSetElementsInternalForces{imat} = find(MESH.MaterialType ==DATAoffline.INITIAL_SET_FOR_ECM_INDEX_MATERIAL(imat) ) ;
    end
    InitialSetElementsInternalForces = cell2mat(InitialSetElementsInternalForces) ;
    ngaus = DATA.MESH.ngaus_STRESS ;
    DATA_ECM.IND_POINTS_CANDIDATES = small2large(InitialSetElementsInternalForces,ngaus) ;
else
    % 11-May-2023
    
    DATAoffline = DefaultField(DATAoffline,'ListElementsInitialCandidates_DECM',[]) ;
    if ~isempty(DATAoffline.ListElementsInitialCandidates_DECM)
        DATA_ECM.IND_POINTS_CANDIDATES = small2large(DATAoffline.ListElementsInitialCandidates_DECM, DATA.MESH.ngaus_STRESS  ) ;
    end
end


DATAoffline.NumberOfCECMpoints = DefaultField(DATAoffline.NumberOfCECMpoints,'INTERNAL_FORCES',[]) ;


DATA_ECM.UseDECMpoints = DATAoffline.UseDECMpoints.INTERNAL_FORCES ;
DATA_ECM.NumberOfCECMpoints = DATAoffline.NumberOfCECMpoints.INTERNAL_FORCES ;

NAMEW_TO_STORE_INTERNAL_FORCEs = [cd,filesep,'MODES',filesep,DATA.NAME_BASE,'CECMpoints',',mat'] ;
DATAoffline = DefaultField(DATAoffline,'LOAD_FROM_MEMORY_ECM_POINTS_INTERNAL_FORCES',0) ;
DATAoffline = DefaultField(DATAoffline,'BINARYFILE_LOAD_FROM_MEMORY_ECM_POINTS_INTERNAL_FORCES',NAMEW_TO_STORE_INTERNAL_FORCEs) ;

if   DATAoffline.LOAD_FROM_MEMORY_ECM_POINTS_INTERNAL_FORCES == 0
    CECM_INTERNAL_FORCES=  CECM_internal_forces(OPERFE,PhiDEF,BasisPone,DATA,MESH,DATA_ECM) ;
    save(DATAoffline.BINARYFILE_LOAD_FROM_MEMORY_ECM_POINTS_INTERNAL_FORCES,'CECM_INTERNAL_FORCES')
else
    disp(['Loading info integration CECM internal forces'])
    load(DATAoffline.BINARYFILE_LOAD_FROM_MEMORY_ECM_POINTS_INTERNAL_FORCES,'CECM_INTERNAL_FORCES')
end

CECM.INTERNAL_FORCES = CECM_INTERNAL_FORCES ;

nmat = length(unique(MESH.MaterialType)) ;
MATLOC = MESH.MaterialType(CECM_INTERNAL_FORCES.setElements);
disp('-----------------------------------------')
disp('Material type for each ECM integration point (internal forces)')
disp('----------------------------------------------------')
for imat = 1:nmat
    nelemLOC  = length(find(MATLOC == imat)) ;
    disp(['Material ',num2str(imat),' = ',num2str(nelemLOC),' points' ])
    
end
disp('-----------------------------------------')


disp('----------------------------------------')
disp('HYPERREDUCTION mass matrix/ BODY FORCES (ASSUMING ONLY GRAVITY IS ACTING)')
disp('------------------------------------------')



DATAoffline = DefaultField(DATAoffline,'errorMASSmatrixDECM',DATAoffline.errorFINT);
DATA_ECM.TOL_SVD_A =DATAoffline.errorMASSmatrixDECM;

DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'MASS_MATRIX',0) ;


DATA_ECM.UseDECMpoints = DATAoffline.UseDECMpoints.MASS_MATRIX ;


DATAoffline = DefaultField(DATAoffline,'LOAD_FROM_MEMORY_ECM_POINTS_MASS_MATRIX',0) ;
DATAoffline = DefaultField(DATAoffline,'BINARYFILE_LOAD_FROM_MEMORY_ECM_POINTS_MASS_MATRIX',NAMEW_TO_STORE_INTERNAL_FORCEs) ;

DATAoffline = DefaultField(DATAoffline,'MassMatrix_ECM_is_dummy_variable',0) ;

 if DATAoffline.MassMatrix_ECM_is_dummy_variable == 1
      
     warning(['Mass matrix and inertial forces not used in this problem'])
     disp(['Accordingly, mass matrix is not accurately computed'])
     DATAoffline.NumberOfCECMpoints.MASS_MATRIX = 2; % DefaultField(DATAoffline.NumberOfCECMpoints,'MASS_MATRIX',2) ;
PhiDEF_mass = PhiDEF(:,1:2) ; 
 else
     DATAoffline.NumberOfCECMpoints = DefaultField(DATAoffline.NumberOfCECMpoints,'MASS_MATRIX',[]) ;
PhiDEF_mass =PhiDEF; 

 end
 DATA_ECM.NumberOfCECMpoints = DATAoffline.NumberOfCECMpoints.MASS_MATRIX ;



if   DATAoffline.LOAD_FROM_MEMORY_ECM_POINTS_MASS_MATRIX == 0
    CECM_MASS_MATRIX=  CECM_mass_matrix(OPERFE,[PhiRB,PhiDEF_mass],MATPRO,DATA,MESH,DATA_ECM) ;
    save(DATAoffline.BINARYFILE_LOAD_FROM_MEMORY_ECM_POINTS_MASS_MATRIX,'CECM_MASS_MATRIX','-append')
else
    disp(['Loading info integration CECM mass matrix'])
    load(DATAoffline.BINARYFILE_LOAD_FROM_MEMORY_ECM_POINTS_MASS_MATRIX,'CECM_MASS_MATRIX')
end

CECM.MASS_MATRIX = CECM_MASS_MATRIX ;



%
%















% diary off

%
%
%
%
% DATAMISC.BASE_FOLDER =   [cd,filesep,'CECMinfo',filesep];
% if  ~exist( DATAMISC.BASE_FOLDER)
%     mkdir( DATAMISC.BASE_FOLDER)
% end
% DATAMISC.ExactIntegral = A'*wFE ;
% DATA_ECM.TOL_SVD_A = DATAoffline.errorFINT ;
% MESH.ngausE = DATA.MESH.ngaus_RHS ;
%
% DATAoffline = DefaultField(DATAoffline,'LoadFromMemory_ECMdata',0) ;
% DATA_ECM = DefaultField(DATA_ECM,'NAMEws_CECMinfo',[cd,'/CECMinfo/','tmp_info.mat']) ;
%
%
% if DATAoffline.LoadFromMemory_ECMdata == 0
%
%     DATA_ECM.ACTIVE =1;  % Enable the continuous empirical cubature method
%     DATA_ECM.TOL_SVD_A =0 ; % Tolerance SVD data matrix  (here the SVD only plays the role of orthogonalization, so we may set it to zero )
%     DATA_ECM.UsePartitionedRandomizedAlgorithm = 0 ;  % Use advanced SVD
%     DATA_ECM.TOL_NewtonRaphson_EliminationPoints = 1e-8 ;  % Tolerance for the Newton-Raphson algorithm
%     DATA_ECM.MaxIterationsNR_ElimPoints = 40 ;  % Number of iteration for each NR step
%     DATA_ECM.SECOND_STAGE_ITERATIONS  = 20 ;   % Number of steps for each weight elimination in the second stage
%
%
%     [CECMoutput,DATA_AUX]= ContinuousECMgen2023(A,xFE,wFE,DATAMISC,MESH,DATA_ECM) ;
%
%     save(DATA_ECM.NAMEws_CECMinfo,'CECMoutput','DATA_AUX') ;
%
%
%
% else
%     load(DATA_ECM.NAMEws_CECMinfo,'CECMoutput','DATA_AUX') ;
%     DATAoffline = DefaultField(DATAoffline,'MAKE_VIDEO_POINTS',0) ;
%     DATA_AUX.AUXVAR = DefaultField(DATA_AUX.AUXVAR,'DATA_from_MAIN',[]) ;
%     DATA_AUX.AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS = DATAoffline.MAKE_VIDEO_POINTS ;
%     DATAoffline = DefaultField(DATAoffline,'BasisIntegrandToPlot_INDEXES',[]) ;
%     DATA_AUX.AUXVAR.DATA_from_MAIN.BasisIntegrandToPlot_INDEXES = DATAoffline.BasisIntegrandToPlot_INDEXES ;
%
% end
%
% VariousPLOTS_CECM(DATA_AUX.DATA_ECM,CECMoutput,DATA_AUX.MESH,DATA_AUX.AUXVAR,...
%     DATA_AUX.VAR_SMOOTH_FE)  ;
%
%
%
% [PHIk_y,~,~]=     EvaluateBasisFunctionALL(CECMoutput.xDECM,DATA_AUX.DATA_ECM,DATA_AUX.VAR_SMOOTH_FE,[]) ;
% b=   DATA_AUX.DATA_ECM.ExactIntegral; % PHI'*W ;
%
%
% bNEW = PHIk_y'*CECMoutput.wDECM ;
% errorINT = norm(bNEW-b)/norm(b)*100;
% disp(['Error DECM rule in integrating basis functions, including interpolation error (%) =',num2str(errorINT)]) ;
%
% %
%
%
% [PHIk_y,~,~]=     EvaluateBasisFunctionALL(CECMoutput.xCECM,DATA_AUX.DATA_ECM,DATA_AUX.VAR_SMOOTH_FE,[]) ;
% b=   DATA_AUX.DATA_ECM.ExactIntegral; % PHI'*W ;
%
%
% bNEW = PHIk_y'*CECMoutput.wCECM ;
% errorINT = norm(bNEW-b)/norm(b)*100;
% disp(['Error CECM rule in integrating basis functions (%) =',num2str(errorINT)]) ;
%
% %
% % %%%
% % ECMdata.COOR = xNEW ;
% % ECMdata.wRED = wNEW ;
% % ECMdata.setPoints = []  ;
% % ECMdata.setElements = ELEMENTS_xNEW ;
% % ECMdata.Ninterpolation = Ninterpolation ;
%
%
%
%
%
%
%
