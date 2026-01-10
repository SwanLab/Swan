function [COORrecons,CNrecons,Materials] = ...
    MeshByRepetitionDOM(COORdom,CNrefNEW,f1NOD,f2NOD,nDOM,MaterialType)


% We have to make several "translated"
% copies of  CNref matrix.
% We begin by renumbering CNref
% if isempty(NODESref)
%     CNrefNEW = CNref ; 
% else
%     
% CNrefNEW = zeros(size(CNref)) ;  % Original connectivity matrix of reference elements 
% for inode = 1:length(NODESref)
%     nodeLOC = NODESref(inode) ;
%     INDnodes = find(CNref==nodeLOC) ;
%     CNrefNEW(INDnodes) = inode ;
% end
% end
% Next we  reorder all matrices by placing as first mesh the one
% corresponding to the reference element.
% REMAINING_meshes = setdiff(1:size(COORnew,2),refMESH) ;
% COORnew = COORnew(:,[refMESH REMAINING_meshes]) ;
% dRECONS =  dRECONS(:,[refMESH REMAINING_meshes]) ;
%-------------------------------------------------------------------
% CONSTRUCTING THE MATRIX OF COORDINATES, CONNECTIVITIES AND MATERIAL
% LIBRARY
% -------------------------------------------------------------------
COORrecons = [COORdom];
CNrecons = [CNrefNEW]  ;
Materials = MaterialType ; 
translationVECTOR = COORdom(f2NOD(1),:)-COORdom(f1NOD(1),:); 
TRANSLATION =0 ; 
nmat = length(unique(Materials)) ; 
for e = 2:nDOM
    % Reference point (local)
    TRANSLATION = TRANSLATION +  translationVECTOR ; 
    COORloc=  COORdom + repmat(TRANSLATION,size(COORdom,1),1) ;
    COORrecons = [COORrecons; COORloc] ;
    CNloc = CNrefNEW +(e-1)*size(COORdom,1) ;
    CNrecons = [CNrecons ; CNloc] ;
          NewMaterial = MaterialType +  (e-1)*nmat ;
    Materials =[Materials; NewMaterial]  ;
    
end