function Nst = AssemblyNGlobalBUB1D(Nelems,MESH,ngaus,DATA,ndim_FINE) 
% Adaption of AssemblyNGlobalBUB.m to 1D problems
% JAHO, 13-mAY-2204, hONEST gREENS, tUSET bARCELONA
% --------------------------------
% AssemblyNGlobalBUB.m is described below ....
% Adaptation to bubble modes (EIFEM) of AssemblyNGlobal.m 
% JAHO, 6-oct-2023, Barcelona, Balmes 185. Spain/13-Oct-2023, Buenas Migas
% (Diag.), Barcelona 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/04_GeneralTheory.mlx
% -------------------------------------------------------------------------
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
% J.A. Hernández, jhortega@cimne.upc.edu  (Initially conceived on 27 Oct. 2015, Barcelona, Spain)
% ------------------------------------------------------------------------------------------------------------------------
%dbstop('14')
if nargin == 0
    load('tmp2.mat')
end

nelem = size(MESH.CN,1) ; 
ndimLEFT = ndim_FINE ; % DOFs per node fine-scale 
ndimRIGHT = MESH.NDOFS_pernode ; % DOFs per node, coarse-scale 
nnode = size(MESH.COOR,1) ; 
nnodeE = size(MESH.CN,2) ; 

m = ndimLEFT*nelem*ngaus ; % Number of rows
n = nnode*ndimRIGHT ;          % Number of columns
nzmaxLOC = m*nnodeE*max(ndimLEFT,ndimRIGHT);   % Maximum number of zeros (number of entries of Belems)
Nst = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst
for inode = 1:nnodeE   % Loop over number of nodes at each element
    setnodes = MESH.CN(:,inode) ;  % Global numbering of the inode-th node of each element
    for idime = 1:ndimRIGHT   % loop over spatial dimensions
        % We employ here the Matlab's built-in "sparse" function --> S = sparse(i,j,s,m,n,nzmax).
        %  This function uses
        % vectors i, j, and s to  generate an m-by-n sparse matrix
        % such that S(i(k),j(k)) = s(k), with space allocated
        % for nzmax nonzeros.
        % --------------------------------------------------------
        % 1) Extracting the column of Nelems corresponding to the idime-th
        % DOF of the inode-th node of each element
        DOFloc = (inode-1)*ndimRIGHT+idime ;  % Corresponding indices
        s = Nelems(:,DOFloc) ;           % Corresponding column
        % -----------------------------------------------------------------
        % 2) Now we determine the position within matrix Bst of each entry
        % of vector "s"
        %  a) Rows
        i = [1:m]' ;  %    (all rows)
        %  b) Columns
        DOFglo = (setnodes-1)*ndimRIGHT+idime ;      %  Global DOFs associated to nodes   "setnodes"
        %  and local DOF "idim". This vector has "nelem" entries,
        % while "i" has
        % "m=nstrain*nelem*ngaus ".
        % Therefore, we have to make
        % nstrain*ngaus tiling
        % copies of DOFglo (see next two lines)
        j = repmat(DOFglo', ndimLEFT*ngaus,1) ; %
        j = reshape(j,length(i),1) ;
        % Assembly process  (as sum of matrix with sparse representation)
        Nst = Nst + sparse(i,j,s,m,n,m) ;
    end
end

Nst = Nst(:,MESH.DOFS_TO_KEEP) ; 
