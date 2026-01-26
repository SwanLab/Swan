function Cdiag = ConvertBlockDiag(C)
% Given the matrices C = [C1;C2..Cn]   where
% size(Ci,1)=size(Ci,2) =m  ConvertBlockDiag returns (in sparse format) the matrix defined
% as 
%  Cdiag = diag(C1,C2,...)
% J.A. Hernández, jhortega@cimne.upc.edu , 27 Oct 2015
%-------------------------------------------------------

if nargin == 0
    C1 = [1 2; 3 4] ; C2 = [5 6; 7 8] ; C3 = [9 10; 15 65];
    C  =[C1;C2;C3] ;  
end
[ii,jj] = IndicesCtang(size(C,1),size(C,2)) ; 
Cdiag =     ConvertCmatSparseMatrix(C,ii,jj)  ;