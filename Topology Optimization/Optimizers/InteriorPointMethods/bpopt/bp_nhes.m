function [h] = bp_nhes(bp,x,s,lam,bL,bU)

n = size(x,2);
m = size(bL,2);

% compute numerical 2nd deriv
dL_base = bp_dL_dx(bp,x,s,lam,bL,bU);
ep = 1e-5;
for i = 1:n,
    xp = x;
    xp(i) = x(i)+ep;
    dL1 = bp_dL_dx(bp,xp,s,lam,bL,bU);
    h(:,i) = [(dL1-dL_base)/ep]';
end
k = 0;
for i = 1:m,
    if (bU(i)>bL(i)),
       k = k + 1;
       sp = s;
       sp(i) = s(i)+ep;
       dL1 = bp_dL_dx(bp,x,sp,lam,bL,bU);
       h(:,n+k) = [(dL1-dL_base)/ep]';
    end
end
