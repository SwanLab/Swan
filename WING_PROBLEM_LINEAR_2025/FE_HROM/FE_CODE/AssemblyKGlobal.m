function K = AssemblyKGlobal(Kelem,nelem,nnodeE,ndim,CN,nnode)
%
% % Given the matrix comprising all element stiffness matrices (Kelem),
% AssemblyKGlobal returns the assembled global stiffness matrix K
% J.A. Hernández, jhortega@cimne.upc.edu  (28 Oct. 2015, Barcelona, Spain)
%dbstop('7')
if nargin == 0
    load('tmp2.mat')
end
m = nnode*ndim ; % Number of rows
n = m ;          % Number of columns
nzmaxLOC = size(Kelem,1)*size(Kelem,2) ;   % Maximum number of zeros (number of entries of Belems)
K = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst

for anod=1:nnodeE % Loop over element nodes (rows)
    a = Nod2DOFelem(anod,ndim,nnodeE,nelem) ;  % ROWS number (in Kelem) for node   "anod" (for all elements)
    for bnod= 1:nnodeE  % Loop over element nodes (columns)
        b = Nod2DOF(bnod,ndim) ; 
        Anod = CN(:,anod) ;  A = Nod2DOF(Anod,ndim) ;  % DOFs in the global K matrix
        Bnod = CN(:,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
        %%%%%
        %  K(A,B) = K(A,B) + Kelem(a,b) ;
        %%%%%
        s = Kelem(a,b) ;  
        s=s(:) ;  nzmax = length(s);
        % Indices "i" and "j"
        i = repmat(A,ndim,1); 
        j = ObtainIndexJassemb(B,ndim) ; 
        %%%%
        
        K = K + sparse(i,j,s,m,n,length(s)) ;
    end
end


% for inode = 1:nnodeE   % Loop over number of nodes at each element
%     setnodes = CN(:,inode) ;  % Global numbering of the inode-th node of each element
%     for idime = 1:ndim   % loop over spatial dimensions
%  
%         % 1) Extracting the column of Belems corresponding to the idime-th
%         % DOF of the inode-th node of each element
%         DOFloc = (inode-1)*ndim+idime ;  % Corresponding indices
%         s = Kelem(:,DOFloc) ;           % Corresponding column
%         % -----------------------------------------------------------------
%         % 2) Now we determine the position within matrix Bst of each entry
%         % of vector "s"
%         %  a) Rows
%         i = [1:length(s)]' ;  %    (all rows)
%         %  b) Columns
%         DOFglo = (setnodes-1)*ndim+idime ;      %  Global DOFs associated to nodes   "setnodes"
%         %  and local DOF "idim". This vector has "nelem" entries,
%         % while "i" has
%         % "m=nstrain*nelem*ngaus ".
%         % Therefore, we have to make
%         % nstrain*ngaus tiling
%         % copies of DOFglo (see next two lines)
%         j = repmat(DOFglo', ndime,1) ; %
%         j = reshape(j,length(i),1) ;
%         % Assembly process  (as sum of matrix with sparse representation)
%         K = K + sparse(i,j,s,m,n,m) ;
%     end
% end