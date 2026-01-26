function S =     ConvertCmatSparseMatrix(C,i,j)  ;
 %S = sparse(i,j,s,m,n,nzmax)  ses
% vectors i, j, and s to
% generate an m-by-n sparse matrix
% such that S(i(k),j(k)) = s(k), with space allocated
% for nzmax nonzeros. Vectors i, j,
% and s are all the same length. Any elements of s that
% are zero are ignored, along with the corresponding values of i and j.
% Any elements of s that have duplicate values of i and j are
% added together.
p = size(C,2);
m = size(C,1) ;
n = m ;
nzmax = m*p ;
s =reshape(C',nzmax,1);
%
S = sparse(i,j,s,m,n,nzmax)  ;