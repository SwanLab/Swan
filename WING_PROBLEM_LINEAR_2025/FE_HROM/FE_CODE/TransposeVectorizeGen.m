function Rt = TransposeVectorizeGen(R,nelem) 

if nargin == 0
    nelem = 2 ; 
    R1 = [1 2 3; 4  5   6] ; 
    R2 = [7 8 9; 10 11 12] ; 
    R = [R1 ;R2] ;
    load('tmp.mat')
end

ndimy = size(R,2) ; 
ndimx = size(R,1)/nelem ;
Rt = zeros(ndimy*nelem,ndimx) ; 
for idim = 1:ndimx
    iglo = idim:ndimx:ndimx*nelem ; 
    for jdim = 1:ndimy
        jglo = jdim:ndimy:ndimy*nelem ; 
        Rt(jglo,idim) = R(iglo,jdim) ; 
    end
end

end