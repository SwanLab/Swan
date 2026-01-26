function [ReactionsFINAL,DISP3D_ROM ]= Training_and_ROM_2Drve(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,MFILENAME,...
    RUN_FE,READ_MESHES,RUN_MODES,RUN_ROM,COMPUTE_MODES_AGAIN,...
    NAME_LOAD_DATA,PLOT_3D_STRUCTURE_CHECK,DATARUN)

% Copy of Training_and_ROM_straightbeam. Multi-scale problems, 2D
% repetitions
% JAHO, 20-July-2018
ReactionsFINAL = []; 

if exist('DATARUN')==0
    DATARUN = [] ;
end
if exist('PLOT_3D_STRUCTURE_CHECK') == 0
    PLOT_3D_STRUCTURE_CHECK = 1;
end


%% Create joints and check geometry
%-----------------------------------
DATAIN.NAMEWS_JOINTS_MESH = ['MODES',filesep,'MESH_',MFILENAME,'.mat'] ;

if READ_MESHES ==1
    cd(FOLDER)
    eval(NAME_INPUT_DATA)
    eval(NAME_LOAD_DATA)
    cd(EXECUTABLE_FOLDER)
    DATAIN.PLOT_3D_STRUCTURE_CHECK = PLOT_3D_STRUCTURE_CHECK ;
    [MESH2D_complete,MESH3D_complete] =GeometryStructureREP2D(MESH2D,MESH3D,DATAIN) ;
    % [MESH2D_complete,MESH3D_complete] =   GENERATE_MESH_JOINTS(NAME_INPUT_DATA,EXECUTABLE_FOLDER,PLOT_3D_STRUCTURE_CHECK) ;
    
    cd(FOLDER)
    if ~exist([FOLDER,filesep,'MODES'],'dir')
        mkdir([FOLDER,filesep,'MODES']) ;
    end
    save(DATAIN.NAMEWS_JOINTS_MESH,'MESH2D_complete','MESH3D_complete') ;
else
    cd(FOLDER)
    load(DATAIN.NAMEWS_JOINTS_MESH,'MESH2D_complete','MESH3D_complete')
end
DATARUN.ndim = size(MESH3D_complete.RVES(1).DATA3D.COOR,2);
DATAIN =DefaultField(DATAIN,'angDOM',[]) ;
% ----------------
% Training SLICES
% -----------------
DATAIN = DefaultField(DATAIN,'TOLERANCE_PERIODIC_NODES_GIVEN',0.1) ; 
cd(FOLDER)

if ~exist('DATA')  % Added 5-April-2022
    DATA = [] ; 
end

if RUN_FE.RVES ==1
    run(NAME_INPUT_DATA) ;
    nrves = length(MESH3D.RVES) ;
    
    for  irve = 1:nrves
        if  isempty(MESH3D_complete.RVES(irve).DATA3D.NODES_CORNERS )%
            % Cellular structure
            DATA = DefaultField(DATA,'BoundaryConditionType','PRESCRIBED_DISPLACEMENTS') ;
       else
            % Continuum structure (plate)
            DATA = DefaultField(DATA,'BoundaryConditionType','PRESCRIBED_DISP_PLATES') ; % Plate BCs
            
            
        end
        DATA= DefaultField(DATA,'FACES_WITH_NON_ZERO_DISPLACEMENTS',[2,3,4]) ; %
        DATARUN.FACES_WITH_NON_ZERO_DISPLACEMENTS = DATA.FACES_WITH_NON_ZERO_DISPLACEMENTS ;
        DATARUN.BoundaryConditionType = DATA.BoundaryConditionType ;
        
        
        
        LOC_SLICE = MESH3D.RVES(irve) ;
        DATARUN.ndim = size(MESH3D_complete.RVES(irve).DATA3D.COOR,2 ) ;
        NameFileMeshLOC = LOC_SLICE.NAME  ;
        LOC_SLICE = DefaultField(LOC_SLICE,'NAME_COARSE',NameFileMeshLOC) ;
        NameFileMeshLOC_coarse = LOC_SLICE.NAME_COARSE ;
        DATAREF = MESH3D.RVES(irve).NAME_DATA_FE ;
        NAME_ROOT_FEloc = [NAME_ROOT_FE.RVES ,num2str(irve),'_'] ;
        DATARUN.angDOM = DATAIN.angDOM;
        DATARUN.TOLERANCE_PERIODIC_NODES_GIVEN =DATAIN.TOLERANCE_PERIODIC_NODES_GIVEN ;
        
        FE_TRAIN_RVE2D(NameFileMeshLOC,NameFileMeshLOC_coarse,...
            DATAREF,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,NUMBER_OF_RVES,PRESCRIBED_VALUE,DATARUN) ;
        
    end
    
    
else
end


cd(FOLDER)
eval(NAME_INPUT_DATA)
cd(EXECUTABLE_FOLDER)
islice = 1;
run(MESH3D.RVES(islice).NAME_DATA_FE) ;
%catch
%  run(NAME_INPUT_DATA)  ;
% run()
% end
FUNinput = DefaultField(FUNinput,'NAME','INPUTS_GENERAL_NONLINEAR') ; 

INPUTS_GENERAL = FUNinput.NAME ;
switch INPUTS_GENERAL
    case 'INPUTS_GENERAL_NONLINEAR'
        DATA.ISNONLINEAR = 1;
        DATARUN.ISNONLINEAR = 1;
    otherwise
        DATA.ISNONLINEAR = 0;
        DATARUN.ISNONLINEAR = 0;
end


 % ---------------------------
% DETERMINATION MODES RVES
% ---------------------------------
if ~exist('TOL_SINGULAR_VALUES_Hqr','var')
    TOL_SINGULAR_VALUES_Hqr.RVES = [] ; 
end

if RUN_MODES.RVES  == 1
    cd(FOLDER)
    run(NAME_INPUT_DATA) ;
    fff =fieldnames(DATAIN);
    for ifield = 1:length(fff)
        DATA.(fff{ifield}) = DATAIN.(fff{ifield})  ;
    end
    
    DATA= DefaultField(DATA,'FACES_WITH_NON_ZERO_DISPLACEMENTS',[1,2,3,4]) ;
    DATARUN.FACES_WITH_NON_ZERO_DISPLACEMENTS = DATA.FACES_WITH_NON_ZERO_DISPLACEMENTS ;
    DATA = DefaultField(DATA,'BoundaryConditionType',[]) ; 
    DATARUN.BoundaryConditionType = DATA.BoundaryConditionType ;
    nrves = length(MESH3D.RVES) ;
    for  irve = 1:nrves
        NAME_MODES  = MESH3D.RVES(irve).NameWSmodes ;
        LOC = MESH3D.RVES(irve) ;
        LOC = DefaultField(LOC,'PROJECTS_TRAIN',[]) ; 
        
        if isempty(LOC.PROJECTS_TRAIN)
            NAME_ROOT_FEloc = [NAME_ROOT_FE.RVES,num2str(irve),'_'] ;            
            LOC = DefaultField(LOC,'NonRVEProjects',[]) ;
            NonRVEProjects = LOC.NonRVEProjects ;
            PROJECTS_TRAIN = [] ; 
        else
            NAME_ROOT_FEloc = [] ; NonRVEProjects = [] ; 
            PROJECTS_TRAIN = LOC.PROJECTS_TRAIN ; 
        end
        
            
        
        GenerateRVEsModes(NAME_MODES,COMPUTE_MODES_AGAIN.RVES,...
            TOL_SINGULAR_VALUES_Hqr.RVES,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,DATA,...
            NonRVEProjects,DATARUN,PROJECTS_TRAIN) ;
    end
end
DISP3D_ROM = [] ; 
if RUN_ROM == 1
    cd(EXECUTABLE_FOLDER)
      [ReactionsFINAL,DISP3D_ROM ]=    ROM_RVE_2Drepetition(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,NAME_LOAD_DATA,MESH2D_complete,...
            MESH3D_complete,DATARUN) ;
  
end
