function  Lbool = AssemblyLbool_gauss(ngaus,CN,nnode,ndim) 

%% GOAL:  Given d (global displacements), and U  (displacements at the Gauss points of  elements CN),
% we wish to contruct two matrices Nall and Lbool such that 
% U = diag(Nall)*Lbool*d  
% Here we focus on Lbool
% All elements are assumed to have the same number of nodes and Gauss
% points
% 
% JAHO, 29-Jun-2021, Cartagena
% --------------------------------------------------------------------------------

%dbstop('15')
if nargin == 0
    load('tmp.mat')
end


nelem = size(CN,1) ;  % Number of elements
nnodeE = size(CN,2) ;  % Number of nodes per element 

% Lbool is a m x n sparse matrix 
m = ndim*nelem*ngaus ; % Number of rows
n = nnode*ndim ;          % Number of columns
nzmaxLOC = m*nnodeE*ndim ;   % Maximum number of zeros (number of entries of Belems)
Lbool = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Lbool

 
for inode = 1:nnodeE   % Loop over number of nodes at each element
    setnodes = CN(:,inode) ;  % Global numbering of the inode-th node of each element
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