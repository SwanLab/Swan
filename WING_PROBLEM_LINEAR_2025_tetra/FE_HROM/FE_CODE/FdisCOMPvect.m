function Fdis = FdisCOMPvect(COOR,CNb,TypeElementB,tracB,idim,CONNECTb) 
%%%%%%%% This function returns the external force vector due to distributed
%%%%%%%% tractions along dimension idim (Vectorized version)
% JAHO, 27-Oct-2015
% ----------------------------------------------------------


%dbstop('9')
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1);  % Number of total nodes
nelemB = size(CONNECTb,1);  % Number of boundary elements (total)
nnodeEb = size(CONNECTb,2) ; % Number of nodes per boundary element
ndim = size(COOR,2); % Number of spatial dimensions
nelemNZ = size(CNb,1); % Number of non-zero components 
Fdis_i = zeros(nnode*ndim,1) ; % Nodal forces caused by traction forces  (dimension idim )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vector Tnod
% --------------
% Tnod_CNb = [tracB(1,1); tracB(1,2); 
%          tracB(2,1); tracB(2,2); 
% ....
%        tracB(nelemB,1)] 
Tnod = sparse(nelemB*nnodeEb,1) ; % Vector of nodal distributed loads 

% Numbering of elements CNb within CONNECTb
% -------------------
% Change 19-July-2017. Why ? 
% -------------------
%    [ setelem]= ElemBnd(CONNECTb,CNb(:)); % elements face "iface"

% ------------------------------------------------------------------------

 setelem = SetElementsBoundary(CONNECTb,CNb) ; 




% Next we assign to vector Tnod the input values  tracB 
for inode = 1:nnodeEb
    ROWS =  (setelem-1)*nnodeEb+inode ; 
    Tnod(ROWS) = tracB(:,inode);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Computing shape function matrices for all boundary elements... (Nst)')
[ NelemB wSTb ] = ComputeNelemBoundALL(COOR,CONNECTb,TypeElementB) ; 
% -------------------------------------------------------------------------
% Assembly of matrix NstB. This is a matrix relating the nodal forces at
% the boundary elements with the nodal forces at the  Gauss points of the
% boundary elements 
ngaus = size(NelemB,1)/nelemB ; 
disp('Assembly of NstB...')
NstB = AssemblyNboundRIGHT(NelemB,nelemB,nnodeEb,ndim-1,ngaus,nnode) ;

%%% Left-operator NstBw'*NstB*Tnod
%%%             
% This is a matrix relating the forces/displacements at all nodes of the
% discretization, with the values of the forces at the Gauss points of the
% boundary elements 
disp('Assembly of NstBw...')
NstBw = AssemblyNboundLEFT(NelemB,nelemB,nnodeEb,ngaus,CONNECTb,nnode) ;
% Diagonal matrix with weights
ndimLOC = 1; 
wDIAGb = CompWeightDiag(wSTb,ndimLOC)  ; 
NstBw = wDIAGb*NstBw ; 
% Boundary forces
% -----------------
NT = NstB*Tnod ; 
%
Fdis_i = NstBw'*NT ;

%% Finally 

Fdis = zeros(nnode*ndim,1) ; 
indices = idim:ndim:nnode*ndim ; 
Fdis(indices) = Fdis_i ;

 
 



% 
% %----------------------------
% for e = 1:nelemB
%     CNloc = CNb(e,:) ;     % Nodes of element "e"
%     tracBe = tracB(e,:)' ; tracBe = tracBe(:) ;   % Nodal values of the distributed load at element "e"
%     Xe = COOR(CNloc,:)' ;% Coordinates of the nodes of element "e"
%     ngaus = length(weig) ; nnodeEb = size(Xe,2)  ;
%     ndimB = size(dershapef,1) ;  Fdis_e = zeros(nnodeEb,1) ;
%     for  g = 1:ngaus
%         % Matrix of derivatives for Gauss point "g"
%         BeXi = dershapef(:,:,g) ;
%         % Matrix of shape functions at point "g"
%         Ne = shapef(g,:) ;
%         %%% Change of coordinates
%         Se = ChangeCoordBnd(Xe,ndimB) ;
%         %%%
%         % Jacobian Matrix
%         Je = Se*BeXi' ;
%         % JAcobian
%         detJe = det(Je) ;
%         %
%         Fdis_e = Fdis_e + weig(g)*detJe*(Ne'*Ne)*tracBe ;
%     end
%     for a=1:nnodeEb
%         Anod = CNb(e,a) ;
%         A = (Anod-1)*ndim +idim ;
%         Fdis_i(A) = Fdis_i(A) + Fdis_i_e(a) ;
%     end
% end