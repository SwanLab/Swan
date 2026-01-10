function [NODES_i,NODES_j] = SetNodesPeriodic1D(DIRICHLET,MESH)
  indLOC = 1; 
    iface = DIRICHLET.FacesPeriodicity1D(indLOC) ; 
    COOR_i = MESH.PROPERTIES_FACES{iface}.COORrelA_global ;  
    NODES_i = MESH.NODES_FACES{iface} ;
    
    
    indLOC = 2; 
    iface = DIRICHLET.FacesPeriodicity1D(indLOC) ; 
    COOR_j = MESH.PROPERTIES_FACES{iface}.COORrelA_global ;  
    NODES_j = MESH.NODES_FACES{iface} ;
    Idx_j = knnsearch(COOR_j,COOR_i) ;
    
     
    COOR_j_i = COOR_j(Idx_j,:) ;
    
    
    errorCOOR = norm(COOR_j_i-COOR_j,'fro') ;
    disp('---------------------------------------------')
    disp(['Checking matching interfaces = ',num2str(1),' and ',num2str(2)]) ;
    disp(['ERROR_match = ',num2str(errorCOOR)])
    if errorCOOR > 1e-6
        error('No periodicity in this problem')
    
        
    end

    % NODES_j(Idx_j)
       % Idx_DOFs_j_i = small2large(Idx_j,ndim) ;
        NODES_j = NODES_j(Idx_j)  ; 