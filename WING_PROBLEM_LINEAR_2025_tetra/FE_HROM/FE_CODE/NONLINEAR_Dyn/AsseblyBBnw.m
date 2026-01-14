function BBnw = AsseblyBBnw(BBnwNA,nstrain,nelem,nnodeE,ndime,ngausE,CONNECT,nnode) ;

if nargin == 0
    load('tmp1.mat')
end
m = nstrain*nelem*ngausE ; 
n = nnode*ndime ; 
nzmax = m*nnodeE*ndime ;  
BBnw = sparse([],[],[],m,n,nzmax);
for inode = 1:nnodeE 
   inodesG = CONNECT(:,inode) ; 
   for idime = 1:ndime 
       iCOL_loc = (inode-1)*ndime+idime ; 
       iCOL_glo = (inodesG-1)*ndime+idime ;
       i = [1:m]' ; 
       j = repmat(iCOL_glo', nstrain*ngausE,1) ;
       j = reshape(j,length(i),1) ;
       s = BBnwNA(:,iCOL_loc) ; 
       BBnw = BBnw + sparse(i,j,s,m,n,m) ; 
   end
end