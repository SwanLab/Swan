function [ReactionsFINAL,DISP3D ]= ROM_RVE_2Drepetition(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,NAME_LOAD_DATA,MESH2D,MESH3D,DATARUN)
if nargin == 0
    load('tmp.mat')
end
% --------------------------------------
% Reduced-order model for a structure made up by repeating a 3D RVE along a
% plane
% JAHO, 24-July-2018 (copy of ROM_straightbeam.m)

% TYPE OF STRUCTURAL ENTITY --> MESH2D.TYPE.INDEX_RVE_JOINT,
% MESH2D.TYPE.ISRVE, MESH2D.TYPE.INDEX3D
% ----------------------------------------------------------
NAME_LOAD_DATA = [FOLDER,filesep,NAME_LOAD_DATA]  ; % Input file
run(NAME_LOAD_DATA)
DISP3D = [] ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd(FOLDER)
%run(NAME_INPUT_DATA_MATERIAL) ;

% Retrieving stiffness matrix and other reduced-order variables.
% -------------------------------------------------------------
TYPE_STRUCTURES = unique(MESH2D.MaterialType) ; % Number of different structural entities
NTYPE = length(TYPE_STRUCTURES) ;
DATAROM_glo = cell(1,NTYPE) ;
DATA_REFMESH_glo = cell(1,NTYPE) ;
DATAIN.NAME_WS_MODES=cell(1,NTYPE) ;
for itypeLOC = 1:NTYPE
    itype = TYPE_STRUCTURES(itypeLOC) ;
    SCRIPT_MODES = MESH2D.PROP(itype).NameWSmodes;
    
    DATAIN.NAME_WS_MODES{itypeLOC} = [FOLDER,filesep,'MODES',filesep,'MODES_',SCRIPT_MODES,'.mat'];
    cd(EXECUTABLE_FOLDER)
    % LOAD_FROM_MEMORY =1 ;ReactionsFINAL
    % if LOAD_FROM_MEMORY ==1
    disp(['Structural entity = ',num2str(itypeLOC)])
    disp('REtrieving OFFLINE data (DATAROM) ...')
    tic
    load(DATAIN.NAME_WS_MODES{itypeLOC},'DATAROM','DATA_REFMESH')
    DATAROM_glo{itype} = DATAROM ;
    DATA_REFMESH_glo{itype} = DATA_REFMESH ;
    toc
    disp('DONE')
end
% else
%     % LOAD BASIS MATRICES
%     % -------------------
%     nBASES_RVE = [] ;
%     load(DATAIN.NAME_WS_MODES,'BASES','DATA_REFMESH','nBASES_RVE') ;
%     DATAROM =  BeamStiffnessMatrix(BASES,DATA_REFMESH,DATAIN,nBASES_RVE) ;
% end
%
% --------------------------------------------------------------------
% 1D MESH
% ---------------------------------------------------------------------

% -----------------------------------------------------------
DATAIN.NAME_INPUT_DATA = NAME_INPUT_DATA ;


% Reduced order model
% -----------------------
DATA_REFMESH_glo{1} = DefaultField(DATA_REFMESH_glo{1},'NODES_CORNERS',[]) ;


DATAIN = DefaultField(DATAIN,'TIME_DISCRETIZATION',[]) ;

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
DATARUN = DefaultField(DATARUN,'OnlyShowPlotReactionForcesROM',0) ;
if isempty(DATAIN.TIME_DISCRETIZATION)
    % One single step, before 12-Jun-2019
    if isempty(DATA_REFMESH_glo{1}.NODES_CORNERS)
        % Open RVEs
        if DATARUN.OnlyShowPlotReactionForcesROM ==0
            [DATAIN,NODES_SNAP,DOFsKEEP,NODES_LINES_ROTATIONS,DISP3D] =     ReducedOrderModel_RVE(DATA_REFMESH_glo,DATAROM_glo,MESH2D,DATAIN,FORCES,DISP,MESH3D,FOLDER,DATARUN) ;
         ReactionsFINAL =    ReactionsPlotROM(DATAIN,NODES_SNAP,MESH2D,DOFsKEEP,NODES_LINES_ROTATIONS) ;
        else
         ReactionsFINAL =    PlotStoredReactions(DATAIN.nameGRAPHS)  ;
        end
    else
        % Close RVEs   (PLATE ELEMENTS )
        ReducedOrderModel_PLATE(DATA_REFMESH_glo,DATAROM_glo,MESH2D,DATAIN,FORCES,DISP,MESH3D,FOLDER,DATARUN) ;
        
    end
    
    
else
    % Several time-steps
    if ~isempty(DATA_REFMESH_glo{1}.NODES_CORNERS)
        errror('Option not implemented yet')
    end
    
    if DATARUN.OnlyShowPlotReactionForcesROM ==0
        [DATAIN,NODES_SNAP,DOFsKEEP,DATA] = ...
            SolveNonLinearRVE_COARSE(DATA_REFMESH_glo,DATAROM_glo,FOLDER,MESH2D,MESH3D,...
            DATARUN,DATAIN,FORCES,DISP,...
            NAME_INPUT_DATA_MATERIAL,EXECUTABLE_FOLDER,NAME_LOAD_DATA) ;
        if ~isempty(DATAIN)
            ReactionsFINAL =  ReactionsPlotROM(DATAIN,NODES_SNAP,MESH2D,DOFsKEEP) ;
        end
    else
        ReactionsFINAL =   PlotStoredReactions(DATAIN.nameGRAPHS)  ;
    end
    
end







