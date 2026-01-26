function Ahinv = UpdateInverseHermitian(Binv,jrow)
% Let B = ([A a]^T*[A a]), where A is a m-by-n full rank  matrix, and "a" a vector of
% m entries. Suppose that we are given the inverse of B (Binv). Function
% UpdateInverseHermitian(Binv,j) returns the value of inv(A^T A) based on the block
% decomposition of % matrix Binv
%
% If the second argument jrow is passed, then it is assumed that B is given
% by
%
%  B = [A(:,1:jrow-1) a A(:,jrow:end)] ;
%
%  J.A. Hdez, 24 March 2017
if nargin == 0
    m = 10; n =4;
    jrow =1 ;
    A = randn(m,n) ; a = randn(m,1) ;
    Bor = [A(:,1:jrow-1) a A(:,jrow:end)] ;
    B = [Bor'*Bor] ;
    Binv = inv(B);
    AhinvREAL = inv(A'*A) 
    
    
end

if nargin ==1
    jrow = size(Binv,2) ;
end

if jrow == size(Binv,2)
    Ahinv = Binv(1:end-1,1:end-1) -  (Binv(1:end-1,end)*Binv(end,1:end-1))/Binv(end,end) ;
else
    AUX = [Binv(:,1:jrow-1)   Binv(:,jrow+1:end)   Binv(:,jrow)] ;
    AUX = [AUX(1:jrow-1,:) ;  AUX(jrow+1:end,:)  ; AUX(jrow,:)] ;
    
    Ahinv = AUX(1:end-1,1:end-1) -  (AUX(1:end-1,end)*AUX(end,1:end-1))/AUX(end,end) ;
    
end



