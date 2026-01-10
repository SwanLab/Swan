function [DATAIN,NODES_SNAP,DOFsKEEP,DATA] = ...
    SolveNonLinearRVE_COARSE(DATA_REFMESH_glo,DATAROM_glo,FOLDER,MESH2D,MESH3D,DATARUN,DATAIN,...
    FORCES,DISP,NAME_INPUT_DATA_MATERIAL,EXECUTABLE_FOLDER,NAME_LOAD_DATA)
if nargin == 0
    load('tmp.mat')
end


% MATERIALS
cd(FOLDER)
run(NAME_INPUT_DATA_MATERIAL)
cd(EXECUTABLE_FOLDER)
if ~exist('J2initialVARIABLES','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/NONLINEAR_Dyn/']);  end


% Non-linear counterpart of ReducedOrderModel_RVE.m  12-Jun-2019
%
[ndimINTF,DATAIN ]= InformationReducedOrderModel(DATAROM_glo,MESH2D,DATAIN) ;
% STEP 1
% ------
% Assembly stiffness matrix
% disp('-------------------------------------')
% disp(['Assembly reduced-order stiffness matrix...'])
% tic
%
% [K,DOFsKEEP ]= AssemblyKrve(DATAROM_glo,MESH2D,DATAIN,ndimINTF) ;
% ndof = size(K,1) ;
%
% disp(['...Done'])
% toc


% 24-June-2019
DATAIN = DefaultField(DATAIN,'AssemblyMethodK_BCB',1) ;  % B'*C*B method for assemblying stiffness matrix  % = 1 ;

TIME_COMP.BMATRICES =tic; 

disp('Assembly coarse-scale global  B-matrix...')
[ wSTs, Bst, BstW, wST,  DATAIN,DOFsKEEP,ASSEMBLY_INFO] = ...
    ComputeBmatricesROMrve(DATAROM_glo,DATA_REFMESH_glo,DATAIN,ndimINTF,MESH2D);
DATAIN.DOFsKEEP = DOFsKEEP;


TIME_COMP.BMATRICES=toc(TIME_COMP.BMATRICES); 
 



disp('-------------------------------------')
% -----------------------------------------

% -------------------------------------
% Assembly vector of external forces
% -----------------------------------
% INTERFACE FORCES
disp('-------------------------------------')
disp(['Assembly reduced-order external force vector...'])
TIME_COMP.EXTERNAL_FORCES =tic; 
 

Ftrac= AssemblyPinterfacesRVE(DATAROM_glo,MESH2D,DATAIN,FORCES,DATA_REFMESH_glo,ndimINTF,DOFsKEEP) ;
TIME_COMP.EXTERNAL_FORCES =toc(TIME_COMP.EXTERNAL_FORCES); 

% -----------------------------------------
% Body and traction forces over elements
% -----------------------------------------
[Fb, fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomainsRVE(DATAROM_glo,MESH2D,DATAIN,FORCES,ndimINTF,DOFsKEEP,DATA_REFMESH_glo) ;

%snapnow;


%F = Fdom + P ;
disp(['Done...'])
disp('-------------------------------------')

 
% -----------------------------------------------------

% Dirichlet boundary conditions
% ------------------------------
TIME_COMP.DIRICHLET_COND =tic ;%(TIME_COMP.DIRICHLET_COND); 

DATAIN.ndimSP = size(DATA_REFMESH_glo{1}.COOR,2) ;
[DOFr,dR,DOFm,Gb] = DirichletBNDCondRVE(DATAROM_glo,MESH2D,DISP,DATAIN,ndimINTF,DOFsKEEP,DATA_REFMESH_glo{1}) ;
% DOFl = 1:size(Bst,2) ;
% DOFl(DOFr) = [] ;
TIME_COMP.DIRICHLET_COND =toc(TIME_COMP.DIRICHLET_COND) ;

% MATERIAL PROPERTIES (THIS SHOULD BE PUT IN A MORE GENERAL FASHION.
% NOW IT IS ONLY VALID FOR J2 PLASTICITY)
% ----------------------------
% This is to construct the arrays need for the nonlinear analysis
% ---> This is for just one slice
% --------------------------------------------
TIME_COMP.nonlinearpr =tic  ;

nstrain = DATAROM_glo{1}.HROMVAR.nstrain ; % Number of strain components (in this case = 4)
setPoints =  DATAROM_glo{1}.HROMVAR.setPoints;  % Selected Gauss points (ECM)
setElements=  DATAROM_glo{1}.HROMVAR.setElements_WITH_REPET;  % Element containing the selected Gauss Points (setElements_WITH_REPET added 29-Apr-2024)
nelem = length(setElements) ;
MaterialType = MESH3D.RVES(1).DATA3D.MaterialType(setElements) ; % Material associated to each element
% Now we repeat MaterialType as many times as number of elements of the
% coarse scale
MaterialType = repmat(MaterialType,size(MESH2D.CN,1),1) ;
nelem = length(MaterialType) ;
ngausE = 1;  % One Gauss point --> One element
% For 1 slice
[DATAIN,PROPMAT,densGLO,CelasGLO] = J2initialVARIABLES(nstrain,DATAIN,nelem,MaterialType,MATERIAL,ngausE) ;
% In a nutshell, if there are NPOINTS integration points, and NCOARSE slices, then the problem is equivalent to one in wich
% the number of elements is equal to NPOINTS*NCOARSE

TypeElement = 'COARSE_SCALE' ;
typePROBLEM = 'pstrain' ; Nst = [] ; M = [] ;

disp('******************************************************************************+')

disp('Solving...')
 
[NODES_SNAP,GAUSS_SNAP,OPERfe,NODESV_PROP,GAUSSV_PROP,DATAIN] =...
    SolveNONLINEAR(Fb,Ftrac,dR,DOFr,MESH2D.COOR,MESH2D.CN,MESH2D.TypeElement,PROPMAT,typePROBLEM,...
    Bst,DOFm,Gb,DATAIN,wST,Nst,BstW,M,ndimINTF,ASSEMBLY_INFO) ;


TIME_COMP.nonlinearpr= toc(TIME_COMP.nonlinearpr) ;
TIME_COMP
total_time = 0 ;
fff = fieldnames(TIME_COMP) ;
for iii=1:length(fff)
    total_time = total_time + TIME_COMP.(fff{iii}) ;
end

disp([' TOTAL TIME NONLINEAR ANALYSIS  = ',num2str(total_time)] )

disp('******************************************************************************+')


save(DATAIN.nameGRAPHS,'TIME_COMP')


%
% % SOLVING SYSTEM OF EQUATIONS
% disp('-------------------------------------')
% tic
% DOFl = 1:size(K,1) ;
% disp(['Solving reduced-order equations  (neq = ',num2str(length(DOFl)),')'])
%
% DOFl(DOFr) = [] ;
% a = zeros(size(K,1),1) ;
% a(DOFr) = aR ;
%
% a(DOFl) = K(DOFl,DOFl)\(F(DOFl)-K(DOFl,DOFr)*aR) ;
% toc
% disp(['...Done'])

% 3D - RECONSTRUCTION PROCESS
% --------------------------------

% PRINTING TIMES (AND FOLDER)
% Different names for the same variable.... bad... Fix it- 10-July-2019



% [~,NAME_LOAD_DATAloc,~ ]= fileparts(NAME_LOAD_DATA) ;
%
% NAME_LOAD_DATAloc = [NAME_LOAD_DATAloc,'_',DATAIN.NAME_INPUT_DATA] ;

NAME_LOAD_DATAloc  = DATAIN.NAME_LOAD_DATAloc   ;

FOLDER_PRINT = [FOLDER,filesep,'GIDPOST',filesep,NAME_LOAD_DATAloc,filesep] ;
if exist(FOLDER_PRINT,'dir')
    rmdir(FOLDER_PRINT,'s') ;
end
mkdir(FOLDER_PRINT) ;
FILES_PRINT = '';
ivar = 1;
NSTEP_REAL = size(NODES_SNAP.U,2) ;
TOTAL_STEPS_TO_PRINT = 1:NSTEP_REAL ;
DATA = DefaultField(DATA,'STEPS_TO_PRINT',TOTAL_STEPS_TO_PRINT ) ;
DATA = DefaultField(DATA,'TIME_TO_PRINT',DATA.STEPS_TO_PRINT ) ;

%DATA = DefaultField(DATA,'TIME_TO_PRINT',TOTAL_STEPS_TO_PRINT) ;
stepsSH = intersect(TOTAL_STEPS_TO_PRINT,[DATA.TIME_TO_PRINT NSTEP_REAL]) ;
stepsSH = sort(stepsSH) ;
% -----------------------------------------

DATAIN  = DefaultField(DATAIN,'ReconstructAndPrintResultsInGid',1) ;

if DATAIN.ReconstructAndPrintResultsInGid ==1
    
    for istepLOC = 1:length(stepsSH)
        istep =     stepsSH(istepLOC) ;
        TIMELOC = DATAIN.TIME_DISCRETIZATION(istep) ;
        a = NODES_SNAP.U(:,istep) ; % Coarse-scale displacement
        stressST = GAUSS_SNAP.stressST(:,istep) ;
        
        % Amplitude self-equilbrated reaction modes
        
        % fextBEAMr and rRB should be multiplied by   DATAIN.FACTOR_TIME_BODY_FORCES
        factMUlT =  cell(size(fextRVEr)) ;
        factMUlT(:) = {DATAIN.FACTOR_TIME_BODY_FORCES(istep)  } ;
        
        fextRVEr_loc = cellfun(@times,fextRVEr,factMUlT,'UniformOutput', false) ;
        rRB_loc = cellfun(@times,rRB,factMUlT,'UniformOutput', false) ;
        fextDOMred_loc = cellfun(@times,fextDOMred,factMUlT,'UniformOutput', false) ;
        
        [rDEF]= AmplitudeReactions_RVE(DATAROM_glo,MESH2D,a,fextRVEr,ndimINTF,DOFsKEEP) ;
        
        
        % Axial, shear forces, torsion and bending moments
        GeneralizedForces= ForcesMoments_RVE(DATAROM_glo,MESH2D,rDEF,rRB,DATA_REFMESH_glo) ;
        % Amplitude displacement modes
        [qDEF,qRB]= AmplitudeDisplacements_RVE(DATAROM_glo,MESH2D,rDEF,fextDOMred,...
            DATA_REFMESH_glo,a,DATAIN,ndimINTF,DOFsKEEP) ;
        % Reconstruction of displacement and stress fields
        % ------------------------------------
        disp('-------------------------------------')
        disp(['Reconstruction of 3D displacement and stresses...'])
        tic
        DATAadd.stressST = stressST;
        [DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA,DATAIN,DATA_REFMESH_glo]= ...
            Displacement_stress_3D_JOINTrve(DATAROM_glo,MESH2D,qDEF,qRB,DATAIN,DATA_REFMESH_glo,DATAadd,DATARUN) ;
        toc
        disp(['...Done'])
        
        
        
        disp('-------------------------------------')
        disp(['Printing results in GID ....'])
        tic
        % Printing post-process file (GID)
        
        % Old impplementation ----just slices. Before July-5th-2018
        % GIDprint_RVE_ROM(MESH2D,DATA_REFMESH,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,STRESS3D,STRESSDATA,...
        %     FORCES_2_PRINT) ;
        
        NameFile_root  = [FOLDER_PRINT,NAME_LOAD_DATAloc,'_',num2str(istep)] ;
        DATAIN.NameFile_msh =  [NameFile_root,'.msh'] ;
        DATAIN.NameFile_res =  [NameFile_root,'.res'] ;
        DATAIN.TIMELOC = TIMELOC ;
        
        
        
        [~,DATAIN] = GIDprint_RVE_ROM(MESH2D,DATA_REFMESH_glo,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,...
            STRESS3D,STRESSDATA,FORCES_2_PRINT,MESH3D,FOLDER,DATARUN,DOFsKEEP,ndimINTF,DATAROM_glo); 
        
        FILES_PRINT = [FILES_PRINT,' ',DATAIN.NameFile_msh,' ', DATAIN.NameFile_res] ;
        
        disp('-------------------------------------')
        toc
        disp(['...Done'])
        
        
        
        
    end
    
    
    DATARUN = DefaultField(DATARUN,'OPEN_GID_ROMRES',0) ;
    clipboard('copy',FOLDER_PRINT)
    disp(FOLDER_PRINT)
    
    %% Plot reaction forces ---along time
    % ReactionsPlotROM(DATAIN,NODES_SNAP,MESH2D,DOFsKEEP) ;
    
    
    
    
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
    
end
