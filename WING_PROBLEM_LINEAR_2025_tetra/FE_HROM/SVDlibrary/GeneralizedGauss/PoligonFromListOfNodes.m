function POLYBND = PoligonFromListOfNodes(LISTofNODES,CONNECTbound)
% LISTofNODES --> List of nodes of a 2D polygon, unsorted
% CONNECTbound --> Connectivity list
% POLYBND --> Sorted lists 
% JAHO, 13-April-2020, 31th day of confinment (COVID-19) 

 LISTofNODES = unique(LISTofNODES) ; 
    % Connectivities 
    POLYBND = zeros(length(LISTofNODES)+1,1) ; 
    POLYBND(1) = LISTofNODES(1) ;  
    nnodeELEM = size(CONNECTbound,2) ;     
    for inode = 1:length(LISTofNODES)        
       [elemBNDloc,inodeELEM] = find(CONNECTbound == POLYBND(inode)) ; 
       % elemBNDloc and inodeELEM should have two elements each 
       % The next point in the polygon is that which is not yet in POLYBND
       ielemLOC = 1; 
       inodeELEMcompl = setdiff(1:nnodeELEM,inodeELEM(ielemLOC)) ; 
       candidateNODE = CONNECTbound(elemBNDloc(ielemLOC),inodeELEMcompl) ; 
       if  ismember(candidateNODE,POLYBND(1:inode)) ;  
           ielemLOC = 2;
           inodeELEMcompl = setdiff(1:nnodeELEM,inodeELEM(ielemLOC)) ;
           candidateNODE = CONNECTbound(elemBNDloc(ielemLOC),inodeELEMcompl) ;
       end
        POLYBND(inode+1) = candidateNODE ;   
     end    
    POLYBND(end) =  POLYBND(1) ;  