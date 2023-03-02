function [dL] = bp_dL_dx(bp,x,s,lam,bL,bU)
n = size(x,2);
m = size(bL,2);

% compute numerical 1st deriv of Lagrangian
L_base = bp_lagrangian(bp,x,s,lam,bL,bU);
ep = 1e-5;
for i = 1:n,
    xp = x;
    xp(i) = x(i)+ep;
    L1 = bp_lagrangian(bp,xp,s,lam,bL,bU);
    dL(i) = [(L1-L_base)/ep];
end
k = 0;
for i = 1:m,
   if (bU(i)>bL(i)),
      k = k + 1;
      sp = s;
      sp(i) = s(i)+ep;
      L1 = bp_lagrangian(bp,x,sp,lam,bL,bU);
      dL(k+n) = [(L1-L_base)/ep];
   end
end

