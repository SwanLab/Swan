function Nst = AssemblyNboundRIGHT(Nelems,nelem,nnodeE,ndim,ngaus,nnode) 
%--------------------------------------------------------------------------
%  AssemblyNboundRIGHT
%
%  Assembles the global **stacked shape function matrix** `Nst` for boundary
%  elements, used to interpolate **nodal tractions** from Gauss-point values
%  on the boundary. This matrix is essential in the weak enforcement of Neumann
%  (natural) boundary conditions in finite element models.
%
%  Mathematically, the matrix relates Gauss-point-level tractions to nodal
%  tractions through a discrete operator:
%
%     t̂(x_gp) ≈ Nst(x_gp) * t̂_node
%
%  which, when integrated and weighted, contributes to the global force vector:
%
%     F̂_trac = Nstᵗ * W * Nst * t̂_node
%
%  INPUTS:
%    - Nelems   : Matrix of evaluated shape functions at Gauss points for each element
%    - nelem    : Number of boundary elements
%    - nnodeE   : Number of nodes per boundary element
%    - ndim     : Spatial dimension (unused here but passed for consistency)
%    - ngaus    : Number of Gauss points per element
%    - nnode    : Total number of nodes in the mesh (unused directly)
%
%  OUTPUT:
%    - Nst : Sparse matrix of size (nelem × ngaus) × (nelem × nnodeE),
%            stacking shape function values for all boundary elements
%
%  ALGORITHM:
%    - Loops over each node position in the local element
%    - For each, assigns Gauss point evaluations into the global sparse matrix
%    - Connectivity is generated virtually (1:nnodeE repeated for each element)
%
%  NOTES:
%    - The matrix returned is purely topological; it maps Gauss-point evaluations
%      to nodal locations for vector-valued boundary quantities.
%    - Used primarily in `NstT_W_N_boundaries` for Neumann condition imposition
%    - Complemented by a left-side operator via `AssemblyNboundLEFT`
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    jhortega@cimne.upc.edu
%    27 Oct. 2015, Barcelona, Spain
%    Comments by ChatGPT4, 13-May-2025
%  Related routines:
%    - `AssemblyNboundLEFT`: builds the corresponding transpose interpolation
%    - `ComputeNelemBoundALL`: evaluates shape functions at boundary Gauss points
%    - Used in computing F̂_trac in Neumann boundary routines
%
%--------------------------------------------------------------------------




% This function returns the global "staked" shape function matrix (boundary
% elements,     the matrix that multiplies tractions forces at nodes )
% J.A. Hernández, jhortega@cimne.upc.edu  (27 Oct. 2015, Barcelona, Spain)
%dbstop('15')
if nargin == 0
    load('tmp1.mat')
end

CONNECT = 1:nnodeE*nelem ;
CONNECT =(reshape(CONNECT', nnodeE,[]))' ; % Local connectibi

m = 1*nelem*ngaus ; % Number of rows
n = nelem*nnodeE ;          % Number of columns
nzmaxLOC = m*nnodeE ;   % Maximum number of zeros
Nst = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst
for inode = 1:nnodeE   % Loop over number of nodes at each element
    setnodes = CONNECT(:,inode) ;  % Global numbering of the inode-th node of each element
    % We employ here the Matlab's built-in "sparse" function --> S = sparse(i,j,s,m,n,nzmax).
    %  This function uses
    % vectors i, j, and s to  generate an m-by-n sparse matrix
    % such that S(i(k),j(k)) = s(k), with space allocated
    % for nzmax nonzeros.
    % --------------------------------------------------------
    % 1) Extracting the column of Nelems corresponding to the idime-th
    % DOF of the inode-th node of each element
    DOFloc = inode ;  % Corresponding indices
    s = Nelems(:,DOFloc) ;           % Corresponding column
    % -----------------------------------------------------------------
    % 2) Now we determine the position within matrix Bst of each entry
    % of vector "s"
    %  a) Rows
    i = [1:m]' ;  %    (all rows)
    %  b) Columns
    DOFglo = setnodes;       %  Global DOFs associated to nodes   "setnodes"
    %  and local DOF "idim". This vector has "nelem" entries,
    % while "i" has
    % "m=nstrain*nelem*ngaus ".
    % Therefore, we have to make
    % nstrain*ngaus tiling
    % copies of DOFglo (see next two lines)
    j = repmat(DOFglo', ngaus,1) ; %
    j = reshape(j,length(i),1) ;
    % Assembly process  (as sum of matrix with sparse representation)
    Nst = Nst + sparse(i,j,s,m,n,m) ;
    
end