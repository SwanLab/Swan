function R = AssemblyRminCOND(CNb,Relem,COOR)
% See ALEX_THESIS_mine.pdf

VECTORIZED = 1 ;

nelem = size(CNb,1) ;
nstrain = size(Relem,1)/nelem ;
ndim = size(COOR,2) ;
nnode =size(COOR,1) ;
nnodeE = size(CNb,2) ;

if VECTORIZED ==1
    m = nstrain ;
    n = ndim*nnode;
    nzmax = prod(size(Relem)) ;
    R = sparse([],[],[],m,n,nzmax) ;
    
    for inode = 1:nnodeE   % Loop over number of nodes at each element
        setnodes = CNb(:,inode) ;  % Global numbering of the inode-th node of each element
        for idime = 1:ndim   % loop over spatial dimensions
            
            % 1) Extracting the column of Relems corresponding to the idime-th
            % DOF of the inode-th node of each element
            DOFloc = (inode-1)*ndim+idime ;  % Corresponding indices
            s = Relem(:,DOFloc) ;           % Corresponding column
            % -----------------------------------------------------------------
            % 2) Now we determine the position within matrix R of each entry
            % of vector "s"
            %  a) Rows
            i = repmat([1:m]',nelem,1) ;  %    (all rows)
            %  b) Columns
            DOFglo = (setnodes-1)*ndim+idime ;      %  Global DOFs associated to nodes   "setnodes"
            %  and local DOF "idim". This vector has "nelem" entries,
            % while "i" has
            % "m=nstrain*nelem*ngaus ".
            % Therefore, we have to make
            % nstrain*ngaus tiling
            % copies of DOFglo (see next two lines)
            j = repmat(DOFglo', nstrain,1) ; %
            j = reshape(j,length(i),1) ;
            % Assembly process  (as sum of matrix with sparse representation)
            R = R + sparse(i,j,s,m,n,length(i)) ;
        end
    end
    
else
    R = zeros(nstrain,ndim*nnode) ;
    
    for ielem = 1:nelem
        CNloc = CNb(ielem,:) ;
        iniE = (ielem-1)*nstrain+1 ;
        finE = ielem*nstrain ;
        Rloc = Relem(iniE:finE,:) ;
        for inode= 1:nnodeE
            for idim = 1:ndim
                DOFloc = (inode-1)*ndim + idim ;
                DOFglo = (CNloc(inode)-1)*ndim + idim ;
                R(:,DOFglo) = R(:,DOFglo) + Rloc(:,DOFloc)  ;
            end
        end
    end
    
end
