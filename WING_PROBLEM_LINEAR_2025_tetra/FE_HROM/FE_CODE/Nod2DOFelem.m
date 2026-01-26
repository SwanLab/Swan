function Aglo = Nod2DOFelem(a,ndim,nnodeE,nelem)
%dbstop('3')
if nargin == 0
    a = 1 ; 
    ndim = 3 ;
    nnodeE = 8 ; 
    nelem = 2 ; 
end


aelem = a:nnodeE:nnodeE*nelem ;  
 
Aglo = Nod2DOF(aelem,ndim) ; 
 