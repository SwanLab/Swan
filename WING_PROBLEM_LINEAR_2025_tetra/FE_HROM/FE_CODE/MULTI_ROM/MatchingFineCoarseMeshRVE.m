function [CONNECTb_faces_coarse,TypeElementB ]= MatchingFineCoarseMeshRVE(NameFileMeshLOC_coarse,COOR)

if nargin == 0
    load('tmp0.mat')
end

    
RVE.NAME =NameFileMeshLOC_coarse ;
DATALOC.InterfacesDoNotMatch =1 ; 
DATA3D = GeometryRVE(RVE,DATALOC) ; % New function, for reading mesh and node faces
NODESfaces = DATA3D.NODES_FACES ;  % All face nodes (labels from GID)

% We only take into account connectivities of the
% faces specified by NODESfaces (by the user, in GID)
% Face 1 and 2 are those normal to the "repetition" direction (plane YZ)
CONNECTb_faces_coarse = cell(1,length(NODESfaces)) ;

TypeElementB = DATA3D.TypeElementB ; 

[zNODES DIST]= knnsearch(COOR,DATA3D.COOR) ;

for iface = 1:length(NODESfaces)
    % Set of nodes 
    % --------------
    NODES_COARSE = NODESfaces{iface} ; % Coarse numbering. Nodes of FACE iface of the coarse mesh 
    [dummy, setBelemLOC]= ElemBnd(DATA3D.CNb,NODES_COARSE); % elements face "iface"
    % Connectivities faces f1 and f2
    CNb_COARSE = DATA3D.CNb(setBelemLOC,:) ;
    
    % REnumber connectivities 
    NODES_FINE = zNODES(NODES_COARSE) ; 
  CONNECTb_faces_coarse{iface} =  RenumberConnectivities(CNb_COARSE,NODES_FINE)  ; 
end


%%%%