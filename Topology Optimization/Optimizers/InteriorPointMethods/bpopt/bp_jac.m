% rearrange Jacobian
function [pd] = bp_jac(bp,x,bL,bU)
pd = bp_jac_stub(bp,x);
m = size(pd,1);
n = size(pd,2);

% add slack variables for inequality constraints
k = 0;
for i = 1:m,
   if(bU(i)>bL(i)),
       k = k + 1;
       pd(i,n+k) = -1;
   end
end
