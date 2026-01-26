function [TypeElementPRINT,MaterialTypegloPRINT,CNprint,COORprint,NAMEMESH] = ...
    PlotGIDFICT_INTERFACES(DATAIN,NAMEMESH,COORprint,DATA_REFMESH,TypeElementPRINT,...
    MaterialTypegloPRINT,CNprint,MESH1D)


if DATAIN.PlotGIDFictitiousInterfaces  == 1
    % Option 22-Dec-2019
    % Plotting in GID interface variables.
    NAMEMESH{end+1} = 'Fict. Interf.' ;
    % MaterialTypegloPRINT, TypeElementPRINT, CNprint, COOR
    % First we have to determine the coordinate of the fictitious interfaces
    nnodesEXIST = size(COORprint,1) ;
    refFACE = 1; itype = 1;
    COORinterface =DATA_REFMESH{itype}.COOR(DATA_REFMESH{itype}.NODES_faces12{refFACE},:) ;
    
    COORinterface  = bsxfun(@minus,COORinterface',DATA_REFMESH{1}.CENTRf1')' ;
    
    CNinterfaceLOC =  RenumberConnectivities( DATA_REFMESH{itype}.CONNECTb{refFACE},1:length(DATA_REFMESH{itype}.NODES_faces12{refFACE})) ;
    TypeElementPRINT{end+1} = DATA_REFMESH{itype}.TypeElementB ;
    
    CNinterface = [] ;
    
    for innode = 1:size(MESH1D.COOR,1) % Loop over 1D nodes
        nSP = size( MESH1D.ROTATIONSint,1) ;
        ifin = nSP*innode;  iini = ifin-nSP+1;
        ROTLOC = MESH1D.ROTATIONSint(:,iini:ifin) ; % Rotation matrix
        if norm(ROTLOC-eye(nSP)) >=1e-10
            error('Option not implemented for curved geometries')
        end
        NEWCOOR = bsxfun(@plus,COORinterface',MESH1D.COOR(innode,:)')' ;
        COORprint = [COORprint; NEWCOOR]  ;
        
        CNinterfaceLOCCC =  CNinterfaceLOC +  nnodesEXIST ;
        CNinterface = [CNinterface ; CNinterfaceLOCCC] ;
        nnodesEXIST = size(COORprint,1) ;
        
        
    end
    MaterialTypegloPRINT{end+1} = 10000*ones(size(CNinterface,1),1) ;
    CNprint{end+1} = CNinterface ;
    
    
    
    
    
    
end