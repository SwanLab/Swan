CNb ={} ; Tnod={} ;
% 1) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 1
%                      %%%%%%%%
idim = 1;
CNb{idim} = [] ;   % In this case, there are no distributed loads in the idim=1 direction
Tnod{idim} = [] ;
% ----
% 2) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 2
%                      %%%%%%%%
idim = 2;
CNb{idim} = [] ;   % In this case, there are no distributed loads in the idim=1 direction
Tnod{idim} = [] ;

% 2) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 3
%                      %%%%%%%%
% Top surface
idim=3 ; 
zmax = max(COOR(BoundaryNodes,3)) ; zmax = zmax(1) ;
nodTOPloc = find(abs(COOR(BoundaryNodes,3)-zmax)<1e-10) ;
NODESb = BoundaryNodes(nodTOPloc) ;
% Connectivity matrix for boundary elements (Neumann boundary)
CNb{idim} = ElemBnd(CONNECTb,NODESb) ;
% Value of the distributed load at the nodes of each boundary element
ty = -20e-3  ; % ;  % MPa
Tnod{idim} = ty*ones(size(CNb{idim})) ;