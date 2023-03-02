function [g] = bp_nobjgrad(bp,x,s)

% compute numerical 1st deriv of the objective function
n = size(x,2);
m = size(s,2);

% compute numerical 1st deriv
obj_base = bp_obj(bp,x);
ep = 1e-5;
for i = 1:n,
    xp = x;
    xp(i) = x(i)+ep;
    obj = bp_obj(bp,xp);
    g(i) = [(obj-obj_base)/ep]';
end
if (m>=1),
    g(n+1:n+m) = 0;
end
