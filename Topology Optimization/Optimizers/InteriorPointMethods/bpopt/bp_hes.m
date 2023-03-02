%% BPOPT Solver: Obtain Hessian (2nd derivatives) of Lagrangian
function [h] = bp_hes(bp,x,s,lam)
n = size(x,2);
ns = size(s,2);
switch bp.prob
    case(0)
        y = apm_details(bp.server,bp.app,x);
        h = y.hes_obj;
        for i = 1:y.neqn,
            h = h + lam(i) * y.hes_eqn{i};
        end
    case(1)
        h = zeros(4,4);
        objfact = 1.0;
        h(1,1) = objfact*2*x(4) + lam(2)*2;                                          %1
        h(2,1) = objfact*x(4) + lam(1)*x(3)*x(4); h(1,2) = h(2,1);                   %2
        h(2,2) = lam(2)*2;                                                           %3
        h(3,1) = objfact*x(4) + lam(1)*x(2)*x(4); h(1,3) = h(3,1);                   %4
        h(3,2) = lam(1)*x(1)*x(4); h(2,3) = h(3,2);                                  %5
        h(3,3) = lam(2)*2;                                                           %6
        h(4,1) = objfact*(2*x(1) + x(2) + x(3)) + lam(1)*x(2)*x(3); h(1,4) = h(4,1); %7
        h(4,2) = objfact*x(1) + lam(1)*x(1)*x(3); h(2,4) = h(4,2);                   %8
        h(4,3) = objfact*x(1) + lam(1)*x(1)*x(2); h(3,4) = h(4,3);                   %9
        h(4,4) = lam(2)*2;                                                           %10
    case(2)
        n = size(x,2);
        h = 2*eye(n);
    case(3)
        h = zeros(2,2);
        h(1,1) = 2;
        h(1,2) = -2;
        h(2,1) = -2;
        h(2,2) = 8;
    case(4)
        h = zeros(2,2);
        h(1,1) = 2 * lam(2);
        h(1,2) = 1 + lam(1);
        h(2,1) = 1 + lam(1);
        h(2,2) = 2 * lam(2);        
    case(5)
        h = zeros(2,2);
        h(1,1) = 2;
        h(1,2) = -2;
        h(2,1) = -2;
        h(2,2) = 8;
    case(6)
        h = zeros(2,2);
        h(1,1) = 2;
        h(2,2) = 4;
    case(7)
        h = zeros(2,2);
        h(1,1) = 2;
        h(2,2) = 2;
    case(8)
        h = zeros(2,2);
        h(1,1) = 2 - 2 * lam(1);
    case(9)
        h(1,1) = 2 + 2 * lam(1);
    case(10)
        h = zeros(2,2);
        h(1,1) = 2;
        h(1,2) = -2;
        h(2,1) = -2;
        h(2,2) = 8;
end

% expand hessian for inequality slacks
for i = 1:ns,
    h(n+i,n+i) = 0.0;
end
