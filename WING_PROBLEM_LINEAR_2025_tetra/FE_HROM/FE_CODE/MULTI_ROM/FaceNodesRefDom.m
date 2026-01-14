function bndDOFS = FaceNodesRefDom(DOMAINVAR,ndim)

idom = 1;
jdom = 1 ;
if isempty(DOMAINVAR.NODES_CORNERS)
    % Cellular structures
    % -------------------    
    faceNODES  = DOMAINVAR.NODES_faces(idom,jdom,:) ;
    faceNODES_I  =faceNODES(:) ;
    bndDOFS_I = cell(size(faceNODES_I)) ;
    for iface = 1:length(faceNODES)
        bndDOFS_I{iface} = small2large(faceNODES{iface},ndim) ;
    end    
    % Therefore, the total list of interface DOFs of domain idom is given
    % by
    bndDOFS = cell2mat(bndDOFS_I) ;    
else
    % Corner DOFs  (4 side elements)
    cornerNODES  = DOMAINVAR.NODES_CORNERS(idom,jdom,:) ;
    cornerNODES_I  =cornerNODES(:) ;
    cornerDOFS_I = cell(size(cornerNODES_I)) ;
    for iface = 1:length(cornerNODES_I)
        cornerDOFS_I{iface} = small2large(cornerNODES_I{iface},ndim) ;
    end
    % Side DOFs
    sideNODES  = DOMAINVAR.NODES_SIDES(idom,jdom,:) ;
    sideNODES_I  =sideNODES(:) ;
    sideDOFS_I = cell(size(sideNODES_I)) ;
    for iface = 1:length(sideNODES_I)
        sideDOFS_I{iface} = small2large(sideNODES_I{iface},ndim) ;
    end
    cornerDOFS = cell2mat(cornerDOFS_I) ;
    sideDOFS = cell2mat(sideDOFS_I)     ;
    
    % Therefore, the total list of interface DOFs of domain idom is given
    % by
    bndDOFS = [cornerDOFS; sideDOFS] ;    
end