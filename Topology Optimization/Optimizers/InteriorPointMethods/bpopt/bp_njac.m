function [J] = bp_njac(bp,x,s,bL,bU)

% compute numerical 1st deriv
n = size(x,2);
m = size(bL,2);

% compute numerical 1st deriv
r_base = bp_res(bp,x,s,bL,bU);
ep = 1e-5;
for i = 1:n,
    xp = x;
    xp(i) = x(i)+ep;
    r = bp_res(bp,xp,s,bL,bU);
    J(:,i) = [(r-r_base)/ep]';
end
k = 0;
for i = 1:m,
    if (bU(i)>bL(i)),
       k = k + 1;
       sp = s;
       sp(i) = s(i)+ep;
       r = bp_res(bp,x,sp,bL,bU);
       J(:,n+k) = [(r-r_base)/ep]';
    end
end
