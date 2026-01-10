CNb ={} ; Tnod={} ;
% 1) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 1
%                      %%%%%%%%
idim = 1;
coorLIM = max(COOR(BoundaryNodes,idim)) ; coorLIM = coorLIM(1) ;
nodes = find(abs(COOR(BoundaryNodes,idim)-coorLIM)<1e-10) ;
NODESb = BoundaryNodes(nodes) ;
% Connectivity matrix for boundary elements (Neumann boundary)
CNb{idim} = ElemBnd(CONNECTb,NODESb) ;
% Value of the distributed load at the nodes of each boundary element
tx = INPUTS_LOC.DISTLOAD.tx  ; % ;  % MPa
Tnod{idim} = tx*ones(size(CNb{idim})) ;
% ----
% 2) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 2
%                      %%%%%%%%
idim = 2;
coorLIM = max(COOR(BoundaryNodes,idim)) ; coorLIM = coorLIM(1) ;
nodes = find(abs(COOR(BoundaryNodes,idim)-coorLIM)<1e-10) ;
NODESb = BoundaryNodes(nodes) ;
% Connectivity matrix for boundary elements (Neumann boundary)
CNb{idim} = ElemBnd(CONNECTb,NODESb) ;
% Value of the distributed load at the nodes of each boundary element
ty = INPUTS_LOC.DISTLOAD.ty  ; % ;  % MPa
Tnod{idim} = ty*ones(size(CNb{idim})) ;

% 2) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 3
%                      %%%%%%%%
% Top surface
idim = 3;
coorLIM = max(COOR(BoundaryNodes,idim)) ; coorLIM = coorLIM(1) ;
nodes = find(abs(COOR(BoundaryNodes,idim)-coorLIM)<1e-10) ;
NODESb = BoundaryNodes(nodes) ;
% Connectivity matrix for boundary elements (Neumann boundary)
CNb{idim} = ElemBnd(CONNECTb,NODESb) ;
% Value of the distributed load at the nodes of each boundary element
tz = INPUTS_LOC.DISTLOAD.tz  ; % ;  % MPa
Tnod{idim} = tz*ones(size(CNb{idim})) ;