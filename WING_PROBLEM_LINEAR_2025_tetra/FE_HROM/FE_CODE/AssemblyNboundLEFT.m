function Nst = AssemblyNboundLEFT(Nelems,nelem,nnodeE,ngaus,CONNECT,nnode) 
%--------------------------------------------------------------------------
%  AssemblyNboundLEFT
%
%  Assembles the global **stacked shape function matrix** `Nst` for boundary
%  elements, used on the **left-hand side** of the boundary integral in the 
%  weak form of Neumann boundary conditions.
%
%  The matrix maps the Gauss-point traction values into the global DOF space,
%  and is used to compute contributions such as:
%
%     F̂_trac = Nstᵗ * W * Nst * t̂_node
%
%  or, in weak form assembly:
%
%     ∫_Γ Nᵗ * t̂(x_gp) dΓ ≈ Nst_leftᵗ * W * t_gp
%
%  INPUTS:
%    - Nelems   : Shape function evaluations at Gauss points for each element
%    - nelem    : Number of boundary elements
%    - nnodeE   : Number of nodes per element
%    - ngaus    : Number of Gauss points per element
%    - CONNECT  : Connectivity matrix of boundary elements (global node IDs)
%    - nnode    : Total number of nodes in the mesh
%
%  OUTPUT:
%    - Nst : Sparse matrix of size (nelem × ngaus) × (nnode × ndim),
%            associating traction field evaluations to DOFs
%
%  ALGORITHM:
%    - Loops over nodes and dimensions (only 1D assumed here)
%    - Places the local contributions from Gauss points into the global matrix
%    - Uses `sparse(i,j,s,m,n)` to build efficiently with memory allocation
%
%  NOTES:
%    - Typically used in computing the left-hand side of the boundary assembly
%      for distributed Neumann loads.
%    - This matrix complements `AssemblyNboundRIGHT`, which maps nodal
%      values to Gauss points (interpolation).
%
%  SPECIAL CASE:
%    - The function assumes `ndim = 1` for simplicity.
%      If vector-valued tractions are to be imposed (e.g., 2D or 3D),
%      the code must be generalized accordingly.
%
%  AUTHOR:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    jhortega@cimne.upc.edu
%    27 Oct. 2015, Barcelona, Spain
%    Comments by ChatGPT3, 13-May-2025
%  SEE ALSO:
%    - `AssemblyNboundRIGHT.m`
%    - `ComputeNelemBoundALL.m` for Nelems generation
%
%--------------------------------------------------------------------------




%This function returns the global "stacked" shape function matrix (boundary
% elements,     the matrix that appears on the left )
% J.A. Hernández, jhortega@cimne.upc.edu  (27 Oct. 2015, Barcelona, Spain)
%dbstop('15')
if nargin == 0
    load('tmp1.mat')
end
ndim= 1; 
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