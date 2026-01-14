function Rt = TransposeVectorize(R) 

if nargin == 0
    R1 = [1 2 ; 3 4] ; 
    R2 = [5 6; 7 8] ; 
    R = [R1 ;R2] ;
end

ndim = size(R,2) ; 
nelem = size(R,1)/ndim ; 
 Rt = zeros(size(R)) ; 
for idim = 1:ndim
    iglo = idim:ndim:ndim*nelem ; 
    for jdim = 1:ndim
        jglo = jdim:ndim:ndim*nelem ; 
        Rt(jglo,idim) = R(iglo,jdim) ; 
    end
end