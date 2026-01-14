function Atop = AtopBot_cal(CNb,nodesZMAX,COOR,TypeElementB,NODESzTOP)
% Compute Atop/Btop matrices
% See Plates_comp_homogenization.pdf

% Connectivity matrix for plane z = zmax
% -------------------------------------- 
nodesINCLUDE =   NODESzTOP ; 
[CNbTOP setBelem]= ElemBnd(CNb,nodesINCLUDE) ;  % Optimize this routine
 
% Compute shape functions for all gauss points  of the top surface
[ NelemB wSTb ] = ComputeNelemBoundALL(COOR,CNbTOP,TypeElementB) ;
% Now we multiply such shape functions by the corresponding weights 
NelemB = bsxfun(@times,NelemB,wSTb) ;
%% Assembly
nelemB = size(CNbTOP,1) ; nnodeEb = size(CNbTOP,2) ;  ngaus = size(NelemB,1)/nelemB ; 
nnode = size(COOR,1) ;
NstBw = AssemblyNboundLEFT(NelemB,nelemB,nnodeEb,ngaus,CNbTOP,nnode) ;
% The above matrix has nnode columns. However,  We are just interested in
% the block matrix corresponding to the nodes nodesZMAX. Thus 
NstBw = NstBw(:,nodesZMAX) ; 
% Finally 
Atop  =sum(NstBw,1) ;