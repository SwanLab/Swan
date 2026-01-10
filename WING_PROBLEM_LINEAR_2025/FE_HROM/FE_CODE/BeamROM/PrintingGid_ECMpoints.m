function [DATAIN] = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HYPERREDUCED_VARIABLES) ;

if nargin == 0
    load('tmp2.mat')
  
end

% Printing both the full and reduced mesh
% ---------------------------------------
DATAIN = DefaultField(DATAIN,'PRINT_MESH_WITH_SELECTED_ELEMENTS_ECM',1) ; % Change 13-Dec-2020

DATAIN = DefaultField(DATAIN,'PRINT_MESH_ECM_POINTS_AS_POINTS',0) ; % Change 15-Feb-2020

% The size of the ECM elements are reduced so that, visually, it looks like
% a point. The idea is to create a new set of points/nodes, and a new reduced mesh
DATAIN  = DefaultField(DATAIN,'COORDINATES_POINTS_TO_PRINT',[]) ; %22-Apr-2020. The coordinates of the points are given 


if  DATAIN.PRINT_MESH_WITH_SELECTED_ELEMENTS_ECM == 1
    %    error('TO BE IMPLEMENTEd')
    % if DATAIN.PLOT_MODES == 1
    % idom = 1;
    %  ELEMS =  DATAOUT.DOMAINVAR.ListElements{idom} ;
    CNref = DATA_REFMESH.CN  ;
    COOR = DATA_REFMESH.COOR  ;
    % LEGENDG = [] ;
    %   NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),LEGENDG ];
    DATA= [] ;
    %  BasisUdef  =[] ;
    
    
    figure(288)
    hold on
    ylabel('Weights (%)')
    xlabel('Points ')
    VOL = sum(HYPERREDUCED_VARIABLES.WdomRED);
    
    WEIGHTS = HYPERREDUCED_VARIABLES.WdomRED/VOL*100 ;
    bar(WEIGHTS)
    
    
    % Name of the mesh fileRpo
    [DIRinp,NAMEFILEinput ]=fileparts(DATAIN.NAME_WS_MODES) ;
    DATAIN = DefaultField(DATAIN,'LABEL_NAME_PROJECT','') ; 
    NAMEFILEinput=  ['ECM_ELEMENTS',DATAIN.LABEL_NAME_PROJECT] ;
    
    ROOTFOLDER = [DIRinp,filesep,'GIDPOST',filesep] ;
    disp(['Writing ECM mesh  in folder: ',ROOTFOLDER]);
    if ~exist(ROOTFOLDER)
        mkdir(ROOTFOLDER) ;
    end
    NameFile_msh = [ROOTFOLDER,filesep,NAMEFILEinput,'' ,'.msh'] ;
    NameFile_res = [ROOTFOLDER,filesep,NAMEFILEinput,'' ,'.res'] ;
    
    
    % Writing mesh file
    %MaterialType = ones(size(CN,1),1) ;
    CNred  = CNref(HYPERREDUCED_VARIABLES.setElements,:) ;
    MatRED = DATA_REFMESH.MaterialType(HYPERREDUCED_VARIABLES.setElements)  ;
    
    
%     HYPERREDUCED_VARIABLES = DefaultField(HYPERREDUCED_VARIABLES,'WEIGHTS_GAUSSIAN_REDUCTION',[]); 
%     
%     if ~isempty(HYPERREDUCED_VARIABLES.WEIGHTS_GAUSSIAN_REDUCTION)
%     end

    
    
    
    if DATAIN.PRINT_MESH_ECM_POINTS_AS_POINTS == 1
        % 15-Feb-2020
        % -----------
        
        if isempty(DATAIN.COORDINATES_POINTS_TO_PRINT)
        disp('Printing mesh for ECM points')
        [COOR,CNredASPOINTS]=  PrintMeshECMPointsAsPoints(COOR,CNred,DATAIN,...
            DATA_REFMESH,HYPERREDUCED_VARIABLES.setElements);
        disp('----------------------------------')
        else
             disp('Printing mesh for ECM points')
        [COOR,CNredASPOINTS]=  PrintMeshECMPointsAsPointsGiven(COOR,CNred,DATAIN,...
            DATA_REFMESH,HYPERREDUCED_VARIABLES.setElements);
        disp('----------------------------------')
        end
        
        
        
        COORprint = COOR ;
        CNprint = {CNredASPOINTS,CNref} ;
        NAME_INPUT_DATA =NAMEFILEinput ;
        MaterialTypegloPRINT = {MatRED,DATA_REFMESH.MaterialType} ;
        TypeElementPRINT = {DATA_REFMESH.TypeElement,DATA_REFMESH.TypeElement} ;
        NAMEMESH = {'REDUCED','FULL'} ;
    else
        COORprint = COOR ;
        CNprint = {CNred,CNref} ;
        NAME_INPUT_DATA =NAMEFILEinput ;
        MaterialTypegloPRINT = {MatRED,DATA_REFMESH.MaterialType} ;
        TypeElementPRINT = {DATA_REFMESH.TypeElement,DATA_REFMESH.TypeElement} ;
        NAMEMESH = {'REDUCED','FULL'} ;
        
    end
    
    
    
    
    
    IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COORprint,CNprint,NAME_INPUT_DATA,MaterialTypegloPRINT,...
        TypeElementPRINT,NAMEMESH);
    %    WEIGHTS = HYPERREDUCED_VARIABLES.WdomRED ;
    
    if length(WEIGHTS)  ~= length(MatRED)
        warning('Weights cannot be plotted')
    else
        
        
        DATA.NAMEMESH = NAMEMESH ;
        DATA.WEIGHTS = WEIGHTS;
        posgp = {[],[]};
        DATAINloc = [] ;
        
        [NameFileRES ]= GidPostProcess_multi_ROM_weights(COORprint,CNprint,TypeElementPRINT, ...
            posgp,NameFile_res,MaterialTypegloPRINT,IND_ELEM_MESHES,...
            DATA_REFMESH,DATA,DATAINloc);
        
    end
    
    
    disp('-----------------------------')
    disp(['To see reduced set of integration points, open']);
    disp(NameFile_msh);
    disp(['--------------------- '])
    
    %
    %     GidPostProcess_onlyECMpoints(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdef,DATA_REFMESH.posgp,...
    %         NAME_MODES,DATA,LEGENDG);
    DATAIN = DefaultField(DATAIN,'MSGPRINT',{});
    DATAIN.MSGPRINT{end+1} = '---------------------------------------------------' ;
    DATAIN.MSGPRINT{end+1} = ['To see reduced set of integration points, open'] ;
    DATAIN.MSGPRINT{end+1} =NameFile_msh ;
    
end
