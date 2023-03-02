%% BPOPT Solver: Obtain Jacobian (1st derivatives) of Equations
function [pd] = bp_jac_stub(bp,x)
switch bp.prob
    case(0)
        y = apm_details(bp.server,bp.app,x);
        pd = y.jac;
    case(1)
        pd(1,1) = x(2)*x(3)*x(4);
        pd(1,2) = x(1)*x(3)*x(4);
        pd(1,3) = x(1)*x(2)*x(4);
        pd(1,4) = x(1)*x(2)*x(3);
        pd(2,1) = 2*x(1);
        pd(2,2) = 2*x(2);
        pd(2,3) = 2*x(3);
        pd(2,4) = 2*x(4);
    case(2)
        n = size(x,2);
        pd = ones(1,n);
    case(3)
        n = size(x,2);
        pd(1,1) = 0.1;
        pd(1,2) = -1;
    case(4)
        pd(1,1) = x(2);
        pd(1,2) = x(1);
        pd(2,1) = 2*x(1);
        pd(2,2) = 2*x(2);
    case(5)
        pd(1,1) = 0.1;
        pd(1,2) = -1.0;
        pd(2,1) = -10.0;
        pd(2,2) = 1.0;        
    case(6)
        pd(1,1) = 2;
        pd(1,2) = 1;
        pd(2,1) = 1;
        pd(2,2) = 2;   
    case(7)
        pd(1,1) = 1;
        pd(1,2) = 1;
    case(8)
        pd(1,1) = -2 * x(1);
        pd(1,2) = -1;
    case(9)
        pd(1,1) = 2 * x(1);
    case(10)
        pd(1,1) = 0.1;
        pd(1,2) = -1.0;
end