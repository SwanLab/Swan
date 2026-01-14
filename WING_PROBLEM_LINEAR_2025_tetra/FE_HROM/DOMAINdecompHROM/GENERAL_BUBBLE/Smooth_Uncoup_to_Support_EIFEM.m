function SMOOTH_from_UNCOUP_TO_SUPPORT = Smooth_Uncoup_to_Support_EIFEM(COOR_SUPPORT,IndicesNEW,nnode_UNCOUP)


ndim  = size(COOR_SUPPORT,2) ; 
nnode_SUPPORT = size(COOR_SUPPORT,1) ;
 P = sparse(nnode_SUPPORT,nnode_UNCOUP) ;

IndicesNEW_aug = [IndicesNEW;(nnode_UNCOUP+1)] ;
for inodeSUPP = 1:nnode_SUPPORT
    INDE_LOC = IndicesNEW_aug(inodeSUPP):(IndicesNEW_aug(inodeSUPP+1)-1) ;
    P(inodeSUPP,INDE_LOC) =  1/length(INDE_LOC)  ;
end
SMOOTH_from_UNCOUP_TO_SUPPORT = sparse(ndim*nnode_SUPPORT,ndim*nnode_UNCOUP) ;
for idim = 1:ndim
    SMOOTH_from_UNCOUP_TO_SUPPORT(idim:ndim:end,idim:ndim:end) = P ;
end