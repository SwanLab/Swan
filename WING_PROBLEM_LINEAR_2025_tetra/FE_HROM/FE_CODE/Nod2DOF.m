function Aglo = Nod2DOF(AnodGLO,ndim)
%dbstop('3')
if nargin == 0
    load('tmp.mat')
end
nnodeE =length(AnodGLO) ;
Aglo = zeros(nnodeE*ndim,1) ;
for inode = 1:nnodeE
    Anod = AnodGLO(inode) ;
    A = [(Anod-1)*ndim+1:Anod*ndim] ;
    if size(A,1) == 1
        A = A' ;
    end
    iini = (inode-1)*ndim+1 ; ifin = inode*ndim ; 
    Aglo(iini:ifin) = A ; 
end
 