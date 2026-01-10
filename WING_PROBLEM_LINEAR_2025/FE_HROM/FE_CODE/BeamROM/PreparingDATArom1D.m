function [MESH1D,DATAIN,DATAROM_glo,DATA_REFMESH_glo,FORCES,ndim,DISP,MATERIAL ]= ...
    PreparingDATArom1D(NAME_INPUT_DATA,NAME_LOAD_DATA,FOLDER,MESH1D,EXECUTABLE_FOLDER,DATAINbef)

%
% 1D -DATA (Supporting, skeleton mesh) --> MESH1D
% ----------------------------------------
% TYPE OF STRUCTURAL ENTITY --> MESH1D.TYPE.INDEX_BEAM_JOINT,
% MESH1D.TYPE.ISBEAM, MESH1D.TYPE.INDEX3D
% ----------------------------------------------------------
NAME_LOAD_DATA_inp = NAME_LOAD_DATA;
NAME_LOAD_DATA = [FOLDER,filesep,NAME_LOAD_DATA]  ; % Input file
run(NAME_LOAD_DATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATAIN = DATAINbef; % 12-May-2020--- I don't know whether this new line may introduce some kind of inconsistency

run(NAME_INPUT_DATA_MATERIAL)
DISP = DefaultField(DISP,'LEFT_END',[] );

if ~isempty(DISP.LEFT_END)
    COOR1_x  = MESH1D.COOR(MESH1D.NODES_POINTS{1},1) ;
    COOR2_x  = MESH1D.COOR(MESH1D.NODES_POINTS{2},1) ;
    IND1 = 1;
    IND2 = 2;
    if COOR2_x < COOR1_x
        IND2 = 1;
        IND1 = 2 ;
    end
    MESH1D.LEFT_END_NODE = MESH1D.NODES_POINTS{IND1}; MESH1D.RIGHT_END_NODE = MESH1D.NODES_POINTS{IND2} ;
else
    % Inputs are introduced in a different way  (directly using points 1 and 2)
end



DATAIN.NODES_ENTITIES  = MESH1D.NODES_POINTS;  % 12-mAY-2020
DATAIN.TableDOFSnode = [] ;  


% Retrieving stiffness matrix and other reduced-order variables.
% -------------------------------------------------------------
TYPE_STRUCTURES = unique(MESH1D.MaterialType) ; % Number of different structural entities
NTYPE = length(TYPE_STRUCTURES) ;
DATAROM_glo = cell(1,NTYPE) ;
DATA_REFMESH_glo = cell(1,NTYPE) ;
DATAIN.NAME_WS_MODES=cell(1,NTYPE) ;
for itypeLOC = 1:NTYPE
    itype = TYPE_STRUCTURES(itypeLOC) ;
    SCRIPT_MODES = MESH1D.PROP(itype).NameWSmodes;
    
    DATAIN.NAME_WS_MODES{itypeLOC} = [FOLDER,filesep,'MODES',filesep,'MODES_',SCRIPT_MODES,'.mat'];
    cd(EXECUTABLE_FOLDER)
    % LOAD_FROM_MEMORY =1 ;
    % if LOAD_FROM_MEMORY ==1
    disp(['Structural entity = ',num2str(itypeLOC)])
    disp('REtrieving OFFLINE data (DATAROM) ...')
    tic
    load(DATAIN.NAME_WS_MODES{itypeLOC},'DATAROM','DATA_REFMESH') ;
    DATAROM_glo{itype} = DATAROM ;
    DATA_REFMESH_glo{itype} = DATA_REFMESH ;
    toc
    disp('DONE')
end

% --------------------------------------------------------------------
% 1D MESH
% ---------------------------------------------------------------------

% -----------------------------------------------------------
DATAIN.NAME_INPUT_DATA = NAME_INPUT_DATA ;


%%%% AS A FIRST APPROACH ---TO BE REVISED IN FUTURE REFINEMENTS OF THE METHOD..., 
% IT IS ASSUMED THAT THE NUMBER OF INTERFACE MODES  IS THE SAME FOR ALL
% DOMAINS (10-OCT-2018)
% WE CHECK THIS CONDITION IN THE SEQUEL

% Checking that all interfaces are the same number of modes 
ndim = [] ; 
for ientity = 1:length(DATAROM_glo)
      V = DATAROM_glo{ientity}.BasisInt; 
     if iscell(V) 
         for iface = 1:length(V)
             ndim(end+1) = size(V{iface},2) ; 
         end
     else
         ndim(end+1) = size(V,2) ; 
     end
end

if any(ndim-ndim(1))
    error('All interfaces should have the same number of modes')
else
    ndim = ndim(1) ; 
end


