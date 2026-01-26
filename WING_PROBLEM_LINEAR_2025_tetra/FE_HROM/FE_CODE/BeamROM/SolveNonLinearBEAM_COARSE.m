function  [DATAIN,NODES_SNAP] = ...
    SolveNonLinearBEAM_COARSE(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,NAME_LOAD_DATA,MESH1D,MESH3D,DATARUN,DATAIN)
if nargin == 0
    load('tmp0.mat')
end
 
 % --------------------------------------
% Reduced-order model for a beam-like structure
% made by repeating  slices. Nonlinear implementation
% JAHO, 29-May-2019
% ----------------------------------------------
% Defining search paths
if ~exist('INPUTS_GENERAL','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/']);  end
if ~exist('SVD_dom','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/MULTILEARN/']);  end
if ~exist('GeometryStructure','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/BeamROM/']);  end
% ------------------

% PREPARING DATA FOR NONLINEAR, ROM ANALYSIS
% ------------------------------------------
[MESH1D,DATAIN,DATAROM_glo,DATA_REFMESH_glo,FORCES,ndimINTF,DISP,MATERIAL ]=...
    PreparingDATArom1D(NAME_INPUT_DATA,NAME_LOAD_DATA,FOLDER,MESH1D,EXECUTABLE_FOLDER,DATAIN) ;
DATAIN.ndimINTF = ndimINTF ; 
DATAIN.npointsECM = length( DATAROM_glo{1}.HROMVAR.setPoints ) ; 

disp('Assembly coarse-scale global  B-matrix...')
[ wSTs, Bst, BstW, wST,  DATAIN] = ...
    ComputeBmatricesROM(DATAROM_glo,DATA_REFMESH_glo,DATAIN,ndimINTF,MESH1D);

% % -------------------------------------
% Assembly vector of external forces
% -----------------------------------
% INTERFACE FORCES
disp('-------------------------------------')
disp(['Assembly coarse-scale external force vector...'])
tic
Ftrac = AssemblyPinterfaces(DATAROM_glo,MESH1D,DATAIN,FORCES,DATA_REFMESH_glo,ndimINTF) ;
% Body and traction forces over elements
[Fb, fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomains(DATAROM_glo,MESH1D,DATAIN,FORCES,ndimINTF,DATA_REFMESH_glo) ;

%F = Fdom + P ;
disp(['Done...'])
disp('-------------------------------------')
toc
% -----------------------------------------------------
% Dirichlet boundary conditions
% ------------------------------
[DOFr,DOFl,dR] = DirichletBNDCondBeam(DATAROM_glo,MESH1D,DISP,DATAIN,ndimINTF) ;


% MATERIAL PROPERTIES (THIS SHOULD BE PUT IN A MORE GENERAL FASHION.
% NOW IT IS ONLY VALID FOR J2 PLASTICITY)
% ----------------------------
% This is to construct the arrays need for the nonlinear analysis
% ---> This is for just one slice
% --------------------------------------------
nstrain = DATAROM_glo{1}.HROMVAR.nstrain ; % Number of strain components (in this case = 4)
setPoints =  DATAROM_glo{1}.HROMVAR.setPoints;  % Selected Gauss points (ECM)
setElements=  DATAROM_glo{1}.HROMVAR.setElements;  % Element containing the selected Gauss Points
nelem = length(setElements) ;
MaterialType = MESH3D.SLICES(1).DATA3D.MaterialType(setElements) ; % Material associated to each element
% Now we repeat MaterialType as many times as number of elements of the
% coarse scale
MaterialType = repmat(MaterialType,size(MESH1D.CN,1),1) ;
nelem = length(MaterialType) ;
ngausE = 1;  % One Gauss point --> One element
% For 1 slice
[DATAIN,PROPMAT,densGLO,CelasGLO] = J2initialVARIABLES(nstrain,DATAIN,nelem,MaterialType,MATERIAL,ngausE) ;
% In a nutshell, if there are NPOINTS integration points, and NCOARSE slices, then the problem is equivalent to one in wich
% the number of elements is equal to NPOINTS*NCOARSE

TypeElement = 'COARSE_SCALE' ;
typePROBLEM = 'pstrain' ; Nst = [] ; M = [] ;
Gb = [] ; DOFm = [] ;  % Variables for imposing affine boundary conditions

disp('Solving...')
ASSEMBLY_INFO = [] ; 
DATAIN.IS_MULTISCALE_ROM_ELEMENT  = 1; 

[NODES_SNAP,GAUSS_SNAP,OPERfe,NODESV_PROP,GAUSSV_PROP,DATAIN] =...
    SolveNONLINEAR(Fb,Ftrac,dR,DOFr,MESH1D.COOR,MESH1D.CN,MESH1D.TypeElement,PROPMAT,typePROBLEM,...
    Bst,DOFm,Gb,DATAIN,wST,Nst,BstW,M,ndimINTF,ASSEMBLY_INFO) ;


% PRINTING TIMES (AND FOLDER)
DATAIN = DefaultField(DATAIN,'TIME_TO_PRINT',1:length(DATAIN.TIME_DISCRETIZATION) ) ;
FOLDER_PRINT = [FOLDER,filesep,'GIDPOST',filesep,NAME_LOAD_DATA,filesep] ;
if exist(FOLDER_PRINT,'dir')
    rmdir(FOLDER_PRINT,'s') ;
end
mkdir(FOLDER_PRINT) ;
FILES_PRINT = '';  
ivar = 1;
NSTEP_REAL = size(NODES_SNAP.U,2) ;
TOTAL_STEPS_TO_PRINT = 1:NSTEP_REAL ;
stepsSH = intersect(TOTAL_STEPS_TO_PRINT,[DATAIN.TIME_TO_PRINT NSTEP_REAL]) ;
stepsSH = sort(stepsSH) ;
% ----------------------------------------- 

for istepLOC = 1:length(stepsSH)
    istep =     stepsSH(istepLOC) ;
    TIMELOC = DATAIN.TIME_DISCRETIZATION(istep) ;
    a = NODES_SNAP.U(:,istep) ; % Coarse-scale displacement
    stressST = GAUSS_SNAP.stressST(:,istep) ;
    
    
    %     % DATAIN.FACTOR_TIME_BODY_FORCES     = DATAIN.TIME_DISCRETIZATION ;
    % DATAIN.FACTOR_TIME_TRACTION_FORCES = DATAIN.TIME_DISCRETIZATION ;
    % DATAIN.FACTOR_TIME_DISPLACEMENTS = DATAIN.TIME_DISCRETIZATION ;
    
    % Plotting interface variables
    DATAIN.PlotInterfaceVariables = 0 ; % Disabled ...
    PlotInterfaceVariablesMOD(DATAIN,a,ndimINTF,MESH1D) ;
    
    % 3D - RECONSTRUCTION PROCESS
    % --------------------------------
    % Amplitude self-equilbrated reaction modes
    
    % fextBEAMr and rRB should be multiplied by   DATAIN.FACTOR_TIME_BODY_FORCES
    factMUlT =  cell(size(fextBEAMr)) ;
    factMUlT(:) = {DATAIN.FACTOR_TIME_BODY_FORCES(istep)  } ;
    
    fextBEAMr_loc = cellfun(@times,fextBEAMr,factMUlT,'UniformOutput', false) ;
    rRB_loc = cellfun(@times,rRB,factMUlT,'UniformOutput', false) ;
    fextDOMred_loc = cellfun(@times,fextDOMred,factMUlT,'UniformOutput', false) ;
    
    [rDEF]= AmplitudeReactions_jbeam(DATAROM_glo,MESH1D,a,fextBEAMr_loc,ndimINTF) ;
    
    % Axial, shear forces, torsion and bending moments
    DATAIN.TIME_STEP_LOCAL = TIMELOC ;
    [GeneralizedForces]= ForcesMoments_diagrams(DATAROM_glo,MESH1D,rDEF,rRB_loc,DATA_REFMESH_glo,ndimINTF,DATAIN) ;
    
    % Amplitude displacement modes
    [qDEF,qRB]= AmplitudeDisplacements(DATAROM_glo,MESH1D,rDEF,fextDOMred_loc,DATA_REFMESH_glo,a,DATAIN,ndimINTF) ;
    
    % Reconstruction of displacement and stress fields
    % ------------------------------------
    disp('-------------------------------------')
    disp(['Reconstruction of 3D displacement and stresses...'])
    tic
    DATAadd.stressST = stressST;
    [DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA,DATAIN,DATA_REFMESH_glo,REACTIONS3D]= ...
        Displacement_stress_3D_JOINTslice(DATAROM_glo,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH_glo,DATARUN,rDEF,rRB_loc,DATAadd) ;
    
    toc
    disp(['...Done'])
    
    
    
    disp('-------------------------------------')
    disp(['Printing results in GID ....'])
    tic
    % Printing post-process file (GID)
    
    % Old impplementation ----just slices. Before July-5th-2018
    % GIDprint_BEAM_ROM(MESH1D,DATA_REFMESH,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,STRESS3D,STRESSDATA,...
    %     FORCES_2_PRINT) ;
    DATAIN.NAME_LOAD_DATA = NAME_LOAD_DATA;
    
    %     ROOTFOLDER = [DIRinp,filesep,'GIDPOST',filesep] ;
    %         if ~exist(ROOTFOLDER)
    %             mkdir(ROOTFOLDER) ;
    %         end
    %         NameFile_msh = [ROOTFOLDER,NameFileMeshHERE,NAMEFILEinput,DATA.POST_LabelFileGid ,'.msh'] ;
    %
    %         NameFile_res= [ROOTFOLDER,NameFileMeshHERE,NAMEFILEinput,DATA.POST_LabelFileGid,'.res'] ;
    
    
    
    NameFile_root  = [FOLDER_PRINT,NAME_LOAD_DATA,'_',num2str(istep)] ;
    DATAIN.NameFile_msh =  [NameFile_root,'.msh'] ;
    DATAIN.NameFile_res =  [NameFile_root,'.res'] ;
    DATAIN.TIMELOC = TIMELOC ;
       GIDprint_BEAM_ROM_JOINTslice(MESH1D,DATA_REFMESH_glo,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,...
        STRESS3D,STRESSDATA,FORCES_2_PRINT,MESH3D,FOLDER,DATARUN,REACTIONS3D)   ;
    
    FILES_PRINT = [FILES_PRINT,' ',DATAIN.NameFile_msh,' ', DATAIN.NameFile_res] ; 
    disp('-------------------------------------')
    toc
    disp(['...Done'])
    
    
end

DATARUN = DefaultField(DATARUN,'OPEN_GID_ROMRES',0) ;
clipboard('copy',FOLDER_PRINT)
disp(FOLDER_PRINT)


if DATARUN.OPEN_GID_ROMRES ==1
%    str = {} ;
 %   str{end+1} = 'Mescape Postprocess';
  %  str{end+1} = 'Mescape ';
    LINENEEE = [' Mescape Files ReadMultiple -rebuildIndex:0',' {',FILES_PRINT,'}'] ;
    clipboard('copy',LINENEEE)
    disp('PASTE INTO THE GID THE STRING contained  IN YOUR CLIPBOARD (see below)')
  disp(LINENEEE)
   % FILE_BATCH = [FOLDER_PRINT,'AUX.txt'] ;
   % fid = fopen(FILE_BATCH,'w') ;
    %             % dbstop('180')
   % for i = 1:length(str)
   %     fprintf(fid,'%s\n',str{i});
   % end
   % fclose(fid ) ;
    
%     str_a =['gid -b ','"',FILE_BATCH,'"',' &'];
%     system(str_a) ;
%     
    
end

%Mescape Files ReadMultiple -rebuildIndex:0 {/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/2D_beam/GIDPOST/LOADS_pbendingLIN/LOADS_pbendingLIN_1.msh /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/2D_beam/GIDPOST/LOADS_pbendingLIN/LOADS_pbendingLIN_1.res /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/2D_beam/GIDPOST/LOADS_pbendingLIN/LOADS_pbendingLIN_2.msh /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/2D_beam/GIDPOST/LOADS_pbendingLIN/LOADS_pbendingLIN_2.res /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/2D_beam/GIDPOST/LOADS_pbendingLIN/LOADS_pbendingLIN_3.msh /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/2D_beam/GIDPOST/LOADS_pbendingLIN/LOADS_pbendingLIN_3.res}

%  
%           
%             %str = GidPostDispl(str) ; % Factor scale
%
%           



