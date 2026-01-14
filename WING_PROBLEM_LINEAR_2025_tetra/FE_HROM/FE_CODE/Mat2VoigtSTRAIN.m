function s = Mat2VoigtSTRAIN(s) ;

if prod(size(s)) == 9

s = [s(1,1); s(2,2); s(3,3) ; 2*s(2,3); 2*s(1,3) ; 2*s(1,2)] ; 
else
   s = [s(1,1); s(2,2);2*s(1,2)] ;  
end


% s = [s(1)  0.5*s(6)  0.5*s(5)
%      0.5*s(6)  s(2)  0.5*s(4)
%      0.5*s(5)  0.5*s(4)  s(3)] ; 
% 
