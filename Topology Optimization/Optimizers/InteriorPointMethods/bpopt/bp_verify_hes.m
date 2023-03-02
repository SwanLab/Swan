function [ok_2nd] = bp_verify_hes(bp,x,s,lam,bL,bU)

% compute numerical 2nd deriv
h_num = bp_nhes(bp,x,s,lam,bL,bU);

% compute exact 2nd deriv
h_exact = bp_hes(bp,x,s,lam);

% evaluate deviation
dev_2nd = max(max(abs(h_num-h_exact)));

tol = 1e-3;
if (dev_2nd>tol),
    fprintf(1,'Hessian error %9.2e exceeds tolerance %9.2e\n',dev_2nd,tol) 
    h_num
    h_exact
    ok_2nd = false;
else
    ok_2nd = true;
end
