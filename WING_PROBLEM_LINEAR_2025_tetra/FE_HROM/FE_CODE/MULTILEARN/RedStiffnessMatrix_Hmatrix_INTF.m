function [KdomREDglo,Hdr,Hrd,Hrr,Hcond] = RedStiffnessMatrix_Hmatrix_INTF(BasisUdef,nDOM,COLUMNS_RVE,KdomRED,BasisRdef,...
    DATAINM,BasisUrb,BasisRrb,alphaBC,NODESbound,ndim)

if nargin == 0
    load('tmp2.mat')
end

if length(BasisUdef)>1 
    error('Option not implemented')
end

% Reduced stiffness matrix 
KdomREDglo = DiagonalGlobalMatrixRVEs(nDOM,KdomRED,COLUMNS_RVE) ;
%%%%%%%%%%%%%%%

% Boundary conditions 
ALPHA = zeros(1,2) ;
if alphaBC(1,1) ==1
    ALPHA(1) = 1;
end
if alphaBC(end,3) ==1
    ALPHA(2) = 1 ;
end

NODESfaces = NODESbound.PLANE ;

%
f1 = small2large(NODESfaces{1},ndim) ;
f2 = small2large(NODESfaces{3},ndim) ;

% Matrix Hdd
BasisR = [BasisRdef{1}] ;
BasisU = [BasisUdef{1}] ;
Hdd = Hmatrices(BasisU,BasisR,f1,f2,nDOM,ALPHA) ; 
% Matrix Hrd
BasisR = [BasisRdef{1}] ;
BasisU = [BasisUrb] ;
Hrd = Hmatrices(BasisU,BasisR,f1,f2,nDOM,ALPHA) ; 
% Matrix Hrr
BasisR = [BasisRrb] ;
BasisU = [BasisUrb] ;
Hrr = Hmatrices(BasisU,BasisR,f1,f2,nDOM,ALPHA) ; 
% Matrix Hdr
BasisR = [BasisRrb] ;
BasisU = [BasisUdef{1}] ;
Hdr = Hmatrices(BasisU,BasisR,f1,f2,nDOM,ALPHA) ;

%%% Condensed matrix 
Hcond = Hdd - Hdr*(Hrr\Hrd) ; 

 