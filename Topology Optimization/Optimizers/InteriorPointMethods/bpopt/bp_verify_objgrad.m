function [ok_1st] = bp_verify_objgrad(bp,x,s)

% compute numerical 1st deriv
g_num = bp_nobjgrad(bp,x,s);

% compute exact 1st deriv
g_exact = bp_objgrad(bp,x,s);

% evaluate deviation
dev_objgrad = max(max(abs(g_num-g_exact)));

tol = 1e-5;
if (dev_objgrad>tol),
    fprintf(1,'Jacobian error %9.2e exceeds tolerance %9.2e\n',dev_objgrad,tol) 
    g_num
    g_exact
    ok_1st = false;
else
    ok_1st = true;
end
