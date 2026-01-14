function Nst = AssemblyNGlobal(Nelems,nelem,nnodeE,ndim,ngaus,CONNECT,nnode) 
%--------------------------------------------------------------------------
% FUNCTION: AssemblyNGlobal
%
% PURPOSE:
%   Constructs the global shape function matrix `Nst` used in finite
%   element formulations for problems with vector-valued unknowns (e.g., displacements).
%   The matrix links the global degrees of freedom to the interpolated field
%   values at Gauss points:
%
%       u(x_gp) ≈ Nst * d
%
%   where:
%     - u(x_gp) are the nodal values of a vector field evaluated at all Gauss points,
%     - d is the global vector of nodal unknowns (e.g., displacements).
%
%   This matrix is block-sparse and has dimensions:
%       Nst ∈ ℝ^{(nelem × ngaus × ndim) × (nnode × ndim)}
%
% INPUTS:
%   - Nelems  : (matrix) Local shape function matrices for all elements
%                        (of size [nelem × ngaus × ndim] × [nnodeE × ndim])
%   - nelem   : (int)    Number of elements
%   - nnodeE  : (int)    Number of nodes per element
%   - ndim    : (int)    Number of spatial dimensions (2 or 3)
%   - ngaus   : (int)    Number of Gauss points per element
%   - CONNECT : (matrix) Element connectivity array (nelem × nnodeE)
%   - nnode   : (int)    Total number of nodes in the mesh
%
% OUTPUT:
%   - Nst     : (sparse matrix) Assembled global matrix of shape functions
%               at all Gauss points, linking global nodal DOFs to field values.
%
% ALGORITHM:
%   For each local node and spatial dimension:
%     1. Extract the local shape function vector associated to that node/direction.
%     2. Map the corresponding global DOFs using the connectivity matrix.
%     3. Replicate those DOF indices across all Gauss points of each element.
%     4. Assemble the triplet (i,j,s) for the sparse matrix, where:
%         - i: row indices (for Gauss points),
%         - j: global DOF indices,
%         - s: shape function values.
%
% REMARKS:
%   - Vectorized for performance using MATLAB's `sparse` function.
%   - Used for consistent mass matrix, internal force evaluation, etc.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   jhortega@cimne.upc.edu
%   CIMNE – Universitat Politècnica de Catalunya
%   27-Oct-2015
%   Comments created by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------




% This function returns the global "staked" shape function matrix given the following
% inputs
%  Nelems: Shape function matrices for all elements  (nelem*ngaus*ndim x nnodeE*ndim)
% nelem: Number of elements
% nnodeE: Number of nodes per element
% ndim: Number of spatial dimensions
% ngaus: Number of gauss per element %
% CONNECT: Array of connectivities
% nnode: Number of nodes
% %
% J.A. Hernández, jhortega@cimne.upc.edu  (27 Oct. 2015, Barcelona, Spain)
%dbstop('14')
if nargin == 0
    load('tmp1.mat')
end
m = ndim*nelem*ngaus ; % Number of rows
n = nnode*ndim ;          % Number of columns
nzmaxLOC = m*nnodeE*ndim ;   % Maximum number of zeros (number of entries of Belems)
Nst = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst
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
        % "m=nstrain*nelem*ngaus ".
        % Therefore, we have to make
        % nstrain*ngaus tiling
        % copies of DOFglo (see next two lines)
        j = repmat(DOFglo', ndim*ngaus,1) ; %
        j = reshape(j,length(i),1) ;
        % Assembly process  (as sum of matrix with sparse representation)
        Nst = Nst + sparse(i,j,s,m,n,m) ;
    end
end