function R = MinCond_Rmatrix(CNbALL,NODESbnd,COOR,TypeElementB)
% Compute R matrix
% See ALEX_THESIS_mine.pdf
if nargin ==0
    load('tmp2.mat')
   % CNbALL = CNb;
    % ---------------------
    COOR = zeros(8,2) ;
    COOR(5:end,:) = [-1 -1; 1 -1 ; 1 1 ; -1 1] ;
    CNbALL = [1 2 ; 2 3; 3 4 ; 4 1; 5 6; 6 7; 7 8; 8 5];
    NODESbnd = [5:8]' ;
    TypeElementB = 'Linear' ;
    
end


[CNb setBelem]= ElemBnd(CNbALL,NODESbnd) ;  %

% Compute Relem for all elements
[ Relem ] = ComputeRelemBoundALL(COOR,CNb,TypeElementB) ;
% ASSEMBLY PROCESS
% ----------------
 R = AssemblyRminCOND(CNb,Relem,COOR) ; 
    
    
    
end







% % Now we multiply such shape functions by the corresponding weights
% NelemB = bsxfun(@times,NelemB,wSTb) ;
% %% Assembly
% nelemB = size(CNbTOP,1) ; nnodeEb = size(CNbTOP,2) ;  ngaus = size(NelemB,1)/nelemB ;
% nnode = size(COOR,1) ;
% NstBw = AssemblyNboundLEFT(NelemB,nelemB,nnodeEb,ngaus,CNbTOP,nnode) ;
% % The above matrix has nnode columns. However,  We are just interested in
% % the block matrix corresponding to the nodes nodesZMAX. Thus
% NstBw = NstBw(:,nodesZMAX) ;
% % Finally
% Atop  =sum(NstBw,1) ;