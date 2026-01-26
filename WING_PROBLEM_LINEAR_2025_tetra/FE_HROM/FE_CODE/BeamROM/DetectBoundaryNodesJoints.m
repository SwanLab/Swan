function [BND ] = DetectBoundaryNodesJoints(CNint,MaterialType1D,itype)
% This function automatically identifies the faces shared by a joint-type domain and their neighboors slices 
% BND.nodes --> Nodes (interfaces )
% BND.slices --> Numbering slices 
%
% See BeamROM.pdf
% JAHO, 29-Dec-2017

% 1D elements of type "itype"
ix = find(MaterialType1D == itype) ;
CN_joint = CNint(ix,:) ;
% Detecting   boundary  nodes  (at least 2) --> BNDnodes
CN_joint = CN_joint(:) ;
BND.nodes = [] ;
for inode = 1:length(CN_joint)
    nodeLOC  =CN_joint(inode) ;
    ix =  find(CN_joint == nodeLOC) ;
    if length(ix) ==1
        BND.nodes = [BND.nodes; CN_joint(ix)] ; 
    end
end
 % Finding connected slices (beam )
 BND.slices = [] ; 
 for inode = 1:length(BND.nodes) % loop over boundary nodes
     nodeLOC = BND.nodes(inode) ;
     [ix iy]= find(CNint==nodeLOC);  % Find boundary nodes in global 1D connectivities
     if  MaterialType1D(ix(1)) ~= itype
         LocBndElemnt = ix(1) ;         
     else
         LocBndElemnt  = ix(2) ;
     end
     BND.slices = [BND.slices; LocBndElemnt] ;
 end