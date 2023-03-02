function [ok_1st] = bp_verify_jac(bp,x,s,bL,bU)

% compute numerical 1st deriv
J_num = bp_njac(bp,x,s,bL,bU);

% compute exact 1st deriv
J_exact = bp_jac(bp,x,bL,bU);

% evaluate deviation
dev_1st = max(max(abs(J_num-J_exact)));

tol = 1e-5;
if (dev_1st>tol),
    fprintf(1,'Jacobian error %9.2e exceeds tolerance %9.2e\n',dev_1st,tol) 
    J_num
    J_exact
    ok_1st = false;
else
    ok_1st = true;
end
