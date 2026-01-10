function W = RepeatWEIGHT(w,n) 
if nargin == 0
    w = [1 2 3 4]'; 
    n = 3 ;
end
w = w' ; 
W =repmat(w,n,1) ; 
W = W(:);