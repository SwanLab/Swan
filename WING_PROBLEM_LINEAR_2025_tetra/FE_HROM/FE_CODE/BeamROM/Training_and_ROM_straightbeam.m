function Training_and_ROM_straightbeam(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,MFILENAME,...
    RUN_FE,GENERATE_JOINTS,RUN_MODES,RUN_ROM,COMPUTE_MODES_AGAIN,...
    NAME_LOAD_DATA,PLOT_3D_STRUCTURE_CHECK,DATARUN)


if exist('DATARUN')==0
    DATARUN = [] ;
end
if exist('PLOT_3D_STRUCTURE_CHECK') == 0
    PLOT_3D_STRUCTURE_CHECK = 1;
end


if ~exist('INPUTS_GENERAL','file') ; addpath(genpath([EXECUTABLE_FOLDER,'/FE_CODE/']));  end
% if ~exist('SVD_dom','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/MULTILEARN/']);  end
% if ~exist('GeometryStructure','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/BeamROM/']);  end
% if ~exist('CurvedDomains','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/MULTI_ROM/']);  end

% ----------------------------------
%% Create joints and check geometry
%-----------------------------------
DATAIN.NAMEWS_JOINTS_MESH = ['MODES',filesep,'MESH_',MFILENAME,'.mat'] ;
RUN_FE = DefaultField(RUN_FE,'SCALE_FACTOR',[]) ;
DATAIN.SCALE_FACTOR  = RUN_FE.SCALE_FACTOR ; % For scaling slices in the X,Y  and Z directions
DATAIN.LOCAL_FOLDER = FOLDER;
if GENERATE_JOINTS ==1
    cd(FOLDER)
    
    % Mod. 13-Nov-2020. Multiple projects
    DATARUN = DefaultField(DATARUN,'INPUTS_NAME_GEOMETRY_MATERIAL_DATA',[]) ;
    
    if  isempty(DATARUN.INPUTS_NAME_GEOMETRY_MATERIAL_DATA)
        eval(NAME_INPUT_DATA)
    else
        nameLOC = 'tmpLOC.mat' ;
        save(nameLOC) ;
        nameLOC = feval(NAME_INPUT_DATA,DATARUN.INPUTS_NAME_GEOMETRY_MATERIAL_DATA,nameLOC) ;
        load(nameLOC) ;
    end
    
    
    
    eval(NAME_LOAD_DATA)
    cd(EXECUTABLE_FOLDER)
    DATAIN.PLOT_3D_STRUCTURE_CHECK = PLOT_3D_STRUCTURE_CHECK ;
    [MESH1D_complete,MESH3D_complete] =GeometryStructure(MESH1D,MESH3D,DATAIN) ;
    % [MESH1D_complete,MESH3D_complete] =   GENERATE_MESH_JOINTS(NAME_INPUT_DATA,EXECUTABLE_FOLDER,PLOT_3D_STRUCTURE_CHECK) ;
    cd(FOLDER)
    if ~exist([FOLDER,filesep,'MODES'],'dir')
        mkdir([FOLDER,filesep,'MODES']) ;
    end
    
    
    save(DATAIN.NAMEWS_JOINTS_MESH,'MESH1D_complete','MESH3D_complete') ;
else
    cd(FOLDER)
    load(DATAIN.NAMEWS_JOINTS_MESH,'MESH1D_complete','MESH3D_complete')
end

DATARUN.ndim = size(MESH1D_complete.COOR,2) ;

DATAIN = DefaultField(DATAIN,'angDOM',[]);
DATAIN = DefaultField(DATAIN,'ELEVATION_Z',[]);

% ----------------
% Training SLICES
% -----------------
cd(FOLDER)
if RUN_FE.SLICES.BEAM ==1
    run(NAME_INPUT_DATA) ;
    nslices = length(MESH3D.SLICES) ;
    DATAIN = DefaultField(DATAIN,'COMBINED_METHOD',1) ;
    for  islice = 1:nslices
        LOC_SLICE = MESH3D.SLICES(islice) ;
        NameFileMeshLOC = LOC_SLICE.NAME  ;
        LOC_SLICE = DefaultField(LOC_SLICE,'NAME_COARSE',NameFileMeshLOC) ;
        NameFileMeshLOC_coarse = LOC_SLICE.NAME_COARSE ;
        DATAREF = MESH3D.SLICES(islice).NAME_DATA_FE ;
        NAME_ROOT_FEloc = [NAME_ROOT_FE.BEAM ,num2str(islice),'_'] ;
        DATARUN.angDOM = DATAIN.angDOM;
        DATARUN.ELEVATION_Z = DATAIN.ELEVATION_Z;
        DATARUN.SCALE_FACTOR = DATAIN.SCALE_FACTOR;
        
        
        if  DATAIN.COMBINED_METHOD  == 0
            warning('this option has not been properly updated...')
            VALUE_FORCE = PRESCRIBED_VALUE.FORCE ;
            VALUE_FORCE = PRESCRIBED_VALUE.MOMENT ;
            [NAME_PROJ_LOC ]=   FE_TRAIN_BEAM_ALLMIXED(NameFileMeshLOC,NameFileMeshLOC_coarse,...
                DATAREF,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,NUMBER_OF_SLICES,PRESCRIBED_VALUE) ;
        else
            
            if ~exist('NUMBER_OF_SLICES','var')
                NUMBER_OF_SLICES = [] ;
                
            end
            
            if ~exist('PRESCRIBED_VALUE','var')
                PRESCRIBED_VALUE = [] ;
            end
            
            
            [NAME_PROJ_LOC ]=   FE_TRAIN_BEAM_COMBINED(NameFileMeshLOC,NameFileMeshLOC_coarse,...
                DATAREF,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,NUMBER_OF_SLICES,PRESCRIBED_VALUE,DATARUN,DATAIN) ;
        end
        
        
    end
    
    
else
end
cd(FOLDER)

if  isempty(DATARUN.INPUTS_NAME_GEOMETRY_MATERIAL_DATA)
    eval(NAME_INPUT_DATA)
else
    nameLOC = 'tmpLOC.mat' ;
    save(nameLOC) ;
    nameLOC = feval(NAME_INPUT_DATA,DATARUN.INPUTS_NAME_GEOMETRY_MATERIAL_DATA,nameLOC) ;
    load(nameLOC) ;
end





cd(EXECUTABLE_FOLDER)
islice = 1;
run(MESH3D.SLICES(islice).NAME_DATA_FE) ;
%catch
%  run(NAME_INPUT_DATA)  ;
% run()
% end
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
% DETERMINATION MODES SLICES
% ---------------------------------
if RUN_MODES.SLICES.BEAM  == 1
    cd(FOLDER)
   % run(NAME_INPUT_DATA) ;
    
    if  isempty(DATARUN.INPUTS_NAME_GEOMETRY_MATERIAL_DATA)
    eval(NAME_INPUT_DATA)
else
    nameLOC = 'tmpLOC.mat' ;
    save(nameLOC) ;
    nameLOC = feval(NAME_INPUT_DATA,DATARUN.INPUTS_NAME_GEOMETRY_MATERIAL_DATA,nameLOC) ;
    load(nameLOC) ;
end

    
    
    
    fff =fieldnames(DATAIN);
    for ifield = 1:length(fff)
        DATA.(fff{ifield}) = DATAIN.(fff{ifield})  ;
    end
    
    nslices = length(MESH3D.SLICES) ;
    DATA.Type_Fluctuation_Modes = 'Separated_by_tests' ;  % New variable to separately treat the fluc. mode of each beam test
    
    DATA = DefaultField(DATA,'UnifiedApproachForModes',1) ;
    DATARUN = DefaultField(DATARUN,'TYPE_BOUNDARY_CONDITIONS','') ;
    switch DATARUN.TYPE_BOUNDARY_CONDITIONS
        case  'PERIODIC_WITH_WARPING' ;
            DATA = DefaultField(DATA,'IncludeWarpingModesAsEssentialBeamModes',1) ; %.IncludeWarpingModesAsEssentialBeamModes = 1;
    end
    
    RUN_MODES = DefaultField(RUN_MODES,'NAME_MODES_FILES',[]) ;
    
    for  islice = 1:nslices
        
        %    if  isempty(RUN_MODES.NAME_MODES_FILES)
        NAME_MODES  = MESH3D.SLICES(islice).NameWSmodes ;
        %   else
        %       NAME_MODES = RUN_MODES.NAME_MODES_FILES ;
        %   end
        
        if ~exist('NAME_ROOT_FE','var')
            NAME_ROOT_FE = [] ;
        end
        
        NAME_ROOT_FE = DefaultField(NAME_ROOT_FE,'BEAM',[]) ;
        
        NAME_ROOT_FEloc = [NAME_ROOT_FE.BEAM,num2str(islice),'_'] ;
        
        LOC = MESH3D.SLICES(islice) ;
        LOC = DefaultField(LOC,'NonBeamProjects',[]) ;
        LOC = DefaultField(LOC,'AdditionalProjects',[]) ;
        
        if DATA.UnifiedApproachForModes == 1
            NonBeamProjects = LOC.AdditionalProjects ; ;
        else
            NonBeamProjects = LOC.NonBeamProjects ;
        end
        
        
        LOC = DefaultField(LOC,'nINTFadditional',[]) ;
        
        DATA.nINTFadditional = LOC.nINTFadditional ;
        
        
        % try
        LOC = DefaultField(LOC,'BeamProjects',[]) ;
        BeamProjects = LOC.BeamProjects ;
        
        TOL_SINGULAR_VALUES_Hqr.SLICES.BEAM = 0.1 ; % OBSOLETE PARAMETER ---16-May-2019
        GenerateBeamModes(NAME_MODES,COMPUTE_MODES_AGAIN.SLICES.BEAM,...
            TOL_SINGULAR_VALUES_Hqr.SLICES.BEAM,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,DATA,...
            NonBeamProjects,DATARUN,BeamProjects) ;
    end
end

% ----------------------------------
% Training JOINTS
% ------------------------------------
DATAIN = DefaultField(DATAIN,'ONLY_DIRICHLET_BCS_JOINTS',1) ;
cd(FOLDER)

if any(RUN_FE.JOINTS)
    run(NAME_INPUT_DATA) ;
    njoints = length(MESH3D.JOINTS) ;
    if length(RUN_FE.JOINTS) ~= njoints
        RUN_FE.JOINTS = ones(njoints,1) ;
    end
    for  ijoint = 1:njoints
        
        if RUN_FE.JOINTS(ijoint)==1
            JOINT = MESH3D.JOINTS(ijoint)  ;
            JOINT = DefaultField(JOINT,'FIXED_FACES',[]) ;
            NameFileMeshLOC = JOINT.NAME  ;
            JOINT =DefaultField(JOINT,'NAME_COARSE',NameFileMeshLOC) ;
            NameFileMeshLOC_coarse = JOINT.NAME_COARSE ;
            DATAREF = JOINT.NAME_DATA_FE ;
            
            if iscell(NAME_ROOT_FE.JOINTS) && length(NAME_ROOT_FE.JOINTS) == njoints
                NAME_ROOT_FEloc = [NAME_ROOT_FE.JOINTS{ijoint},'_'] ;
            else
                NAME_ROOT_FEloc = [NAME_ROOT_FE.JOINTS ,num2str(ijoint),'_'] ;
            end
            
            
            if  DATAIN.ONLY_DIRICHLET_BCS_JOINTS  == 0
                
                FE_TRAIN_JOINTS_ALLMIXED(NameFileMeshLOC,NameFileMeshLOC_coarse,DATAREF,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,...
                    MESH3D_complete,MESH1D_complete,ijoint,NUMBER_OF_SLICES_JOINTS,VALUE_FORCE_JOINTS,...
                    VALUE_MOMENT_JOINTS) ;
            else
                
                if isempty(JOINT.FIXED_FACES)
                    FE_TRAIN_JOINTS_ONLYDISP(NameFileMeshLOC,NameFileMeshLOC_coarse,DATAREF,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,...
                        MESH3D_complete,MESH1D_complete,ijoint,NUMBER_OF_SLICES_JOINTS,PRESCRIBED_VALUE) ;
                else
                    FE_TRAIN_JOINTS_OD_FIXED(NameFileMeshLOC,NameFileMeshLOC_coarse,DATAREF,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,...
                        MESH3D_complete,MESH1D_complete,ijoint,NUMBER_OF_SLICES_JOINTS,PRESCRIBED_VALUE,JOINT.FIXED_FACES) ;
                end
                
            end
        end
    end
end


% ---------------------------
% DETERMINATION MODES JOINTS
% ---------------------------------
if any(RUN_MODES.JOINTS)
    
    njoints = length(MESH3D.JOINTS) ;
    if length(RUN_MODES.JOINTS) ~= njoints
        RUN_MODES.JOINTS = ones(njoints,1) ;
        COMPUTE_MODES_AGAIN.JOINTS  =ones(njoints,1) ;
    end
    cd(FOLDER)
    run(NAME_INPUT_DATA) ;
    njoints = length(MESH3D.JOINTS) ;
    DATA.IS_A_JOINT = 1 ;
    DATA.SLICES_SEL= [] ;
    for  ijoint = 1:njoints
        if RUN_MODES.JOINTS(ijoint)==1
            NAME_MODES  = MESH3D.JOINTS(ijoint).NameWSmodes ;
            
            if iscell(NAME_ROOT_FE.JOINTS) && length(NAME_ROOT_FE.JOINTS) == njoints
                NAME_ROOT_FEloc = [NAME_ROOT_FE.JOINTS{ijoint},'_'] ;
            else
                NAME_ROOT_FEloc = [NAME_ROOT_FE.JOINTS ,num2str(ijoint),'_'] ;
            end
            NonBeamProjects = [] ;
            GenerateBeamModes(NAME_MODES,COMPUTE_MODES_AGAIN.JOINTS(ijoint),...
                TOL_SINGULAR_VALUES_Hqr.JOINTS,NAME_ROOT_FEloc,EXECUTABLE_FOLDER,DATA,NonBeamProjects) ;
        end
    end
end

if RUN_ROM == 1
    run(NAME_LOAD_DATA) ;
    DATAIN = DefaultField(DATAIN,'TIME_DISCRETIZATION',[]) ;
    cd(EXECUTABLE_FOLDER)
    
    
    if isempty(DATAIN.TIME_DISCRETIZATION)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Linear elastic range, just one step
        % ------------------------------------- Version before 9-Oct-2018
        ROM_straightbeam(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,NAME_LOAD_DATA,MESH1D_complete,...
            MESH3D_complete,DATARUN) ;
    else
        % ----------------
        % Nonlinear regime
        % ----------------
        DATARUN = DefaultField(DATARUN,'OnlyShowPlotReactionForces',0) ;
        % ----------------------
        % For plotting graphs
        % ---------------------
        FOLDER_PRINT = [FOLDER,filesep,'GRAPHS'] ;
        if ~exist(FOLDER_PRINT,'dir')
            mkdir(FOLDER_PRINT)
        end
        [~,NAME_LOAD_DATAloc,~ ]= fileparts(NAME_LOAD_DATA) ;
        NAME_LOAD_DATAloc = [NAME_LOAD_DATAloc,'_',NAME_INPUT_DATA] ;
        DATAIN.NAME_LOAD_DATAloc = NAME_LOAD_DATAloc ;
        DATAIN.nameGRAPHS = [FOLDER_PRINT,filesep,DATAIN.NAME_LOAD_DATAloc,'.mat'] ;
        
        if DATARUN.OnlyShowPlotReactionForces  ==0
            %  [DATAIN,NODES_SNAP,DOFsKEEP,DATA] = ...
            [DATAIN,NODES_SNAP] =  SolveNonLinearBEAM_COARSE(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,NAME_LOAD_DATA,MESH1D_complete,...
                MESH3D_complete,DATARUN,DATAIN) ;
            %   if ~isempty(DATAIN)
            DOFsKEEP = [] ;
            ReactionsPlotROM1D(DATAIN,NODES_SNAP,MESH1D) ;
            %  end
        else
            PlotStoredReactions(DATAIN.nameGRAPHS)  ;
        end
        
        
        
        
    end
    % -------------------------------------
end
