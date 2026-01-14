function PlotModesFLUCT_interfaceLevel(Vrb,Ufluc,faceDOFS,DATAinputLOC,REFMESH,cnINTF)
if nargin == 0
    load('tmp.mat')
end

FOLDER_MODES =[DATAinputLOC.FOLDER_PROJECT,filesep,'MODES_INTF',filesep] ;
if ~exist(FOLDER_MODES)
    mkdir(FOLDER_MODES) ;
end

for I = 1:length(Vrb)
    disp(['Interface = ',num2str(I),', post-processing disp. modes'])
    
    if  ~isempty(Ufluc)  
        
        MODES=  [ Vrb{I}/norm(Vrb{I},'fro'),Ufluc{I}/norm(Ufluc{I},'fro')] ;
        
    else
        MODES=  [ Vrb{I}/norm(Vrb{I},'fro')] ;
    end
    
    
    % Associated domain
    % -----------------------
     idoms = find(I == cnINTF) ;
    [INDEX_DOM,INDEX_FACE] =   ind2sub(size(cnINTF),idoms) ;
    INDEX_DOM = INDEX_DOM(1) ; INDEX_FACE = INDEX_FACE(1); 
    
    NODES_faces  = REFMESH{INDEX_DOM}.NODES_faces12{INDEX_FACE} ; 
    COOR= REFMESH{INDEX_DOM}.COOR(NODES_faces,:) ;  
    CNbdom = REFMESH{INDEX_DOM}.CONNECTb{INDEX_FACE} ; 
    CN  =  RenumberConnectivities( CNbdom,1:length(NODES_faces) );
     
    
    NAME_MODES = ['Intf. disp. modes'] ;
    DATA = [] ;
    LEGENDG = 'Intf. disp. modes' ;
    MSG = [] ;
    nameLCO =  [FOLDER_MODES,'InterfDispModes','_INTF_',num2str(I)] ;
    DATA.NameFile_msh =[nameLCO,'.msh'];
    DATA.NameFile_res = [nameLCO,'.res'];
    
    NAME_MODES = [] ;
    
    GidPostProcessModes_dom(COOR,CN,REFMESH{INDEX_DOM}.TypeElementB,MODES,[],...
        NAME_MODES,DATA,LEGENDG,MSG);
    
    
    %       if DATAIN.ISGAUSS == 0
    %    MSG  =   GidPostProcessModes_dom(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdef,DATA_REFMESH.posgp,...
    %         NAME_MODES,DATA,LEGENDG,MSG);
    %     else
    %       MSG =    GidPostProcessModes_domGAUSS(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdef,DATA_REFMESH.posgp,...
    %         NAME_MODES,DATA,LEGENDG,MSG);
    %     end
    
end