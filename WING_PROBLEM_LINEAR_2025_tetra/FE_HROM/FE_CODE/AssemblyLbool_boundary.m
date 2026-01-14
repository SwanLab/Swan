function Lbool = AssemblyLbool_boundary(ngaus,CONNECT,nnode,ndim) 

%% GOAL:  Given d (global displacements), and U  (displacements at the Gauss points of boundary elements CONNECT),
% we wish to contruct two matrices NbndALL and Lbool such that 
% U = diagBLOCK(NbndALL)*Lbool*d v
% Here we focus on the coLboolruction of the sparse boolean operator  Lbool 
% All elements are assumed to have the same number of nodes
% 
% JAHO, 29-Jun-2021, Cartagena
% --------------------------------------------------------------------------------

%dbstop('15')
if nargin == 0
    load('tmp1.mat')
end

nelem = size(CONNECT,1) ; 
nnodeE = size(CONNECT,2) ; 

m = ndim*nelem*ngaus ; % Number of rows
n = nnode*ndim ;          % Number of columns
nzmaxLOC = m*nnodeE*ndim ;   % Maximum number of zeros (number of entries of Belems)
Lbool = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst

Nelems 

for inode = 1:nnodeE   % Loop over number of nodes at each element
    setnodes = CONNECT(:,inode) ;  % Global numbering of the inode-th node of each element
    for idime = 1:ndim   % loop over spatial dimensions
        % We employ here the Matlab's built-in "sparse" function --> S = sparse(i,j,s,m,n,nzmax).
        %  This function uses
        % vectors i, j, and s to  generate an m-by-n sparse matrix
        % such that S(i(k),j(k)) = s(k), with space allocated
        % for nzmax nonzeros.
        % --------------------------------------------------------
        % 1) Extracting the column of Nelems corresponding to the idime-th
        % DOF of the inode-th node of each element
        DOFloc = (inode-1)*ndim+idime ;  % Corresponding indices
        s = Nelems(:,DOFloc) ;           % Corresponding column
        % -----------------------------------------------------------------
        % 2) Now we determine the position within matrix Bst of each entry
        % of vector "s"
        %  a) Rows
        i = [1:m]' ;  %    (all rows)
        %  b) Columns
        DOFglo = (setnodes-1)*ndim+idime ;      %  Global DOFs associated to nodes   "setnodes"
        %  and local DOF "idim". This vector has "nelem" entries,
        % while "i" has
        % "m=Lboolrain*nelem*ngaus ".
        % Therefore, we have to make
        % Lboolrain*ngaus tiling
        % copies of DOFglo (see next two lines)
        j = repmat(DOFglo', ndim*ngaus,1) ; %
        j = reshape(j,length(i),1) ;
        % Assembly process  (as sum of matrix with sparse representation)
        Lbool = Lbool + sparse(i,j,s,m,n,m) ;
    end
end