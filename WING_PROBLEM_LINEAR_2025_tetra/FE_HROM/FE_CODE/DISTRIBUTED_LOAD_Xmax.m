CNb =cell(3,1); Tnod=cell(3,1) ;
 
% 
% % 1) List of boundary nodes whose corresponding boundary elements have distributed loads
% % along  direction --> idim = 1
% %                      %%%%%%%%
% idim = 1;
% CNb{idim} = [] ;   % In this case, there are no distributed loads in the idim=1 direction
% Tnod{idim} = [] ;
% % ----
% % 2) List of boundary nodes whose corresponding boundary elements have distributed loads
% % along  direction --> idim = 2
% %                      %%%%%%%%
% idim = 2;
% CNb{idim} = [] ;   % In this case, there are no distributed loads in the idim=1 direction
% Tnod{idim} = [] ;

% 2) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 3
%                      %%%%%%%%
% xmax surface
idim=1 ; 
xmax = max(COOR(BoundaryNodes,1)) ; xmax = xmax(1) ;
nodes = find(abs(COOR(BoundaryNodes,1)-xmax)<1e-10) ;
NODESb = BoundaryNodes(nodes) ;
% Connectivity matrix for boundary elements (Neumann boundary)
CNb{idim} = ElemBnd(CONNECTb,NODESb) ;
% Value of the distributed load at the nodes of each boundary element
ty = 100 ;  % kN/m
Tnod{idim} = ty*ones(size(CNb{idim})) ;