%% BPOPT Solver: Obtain equation residuals
function [c] = bp_res_stub(bp,x)
switch bp.prob
    case(0)
        y = apm_details(bp.server,bp.app,x);
        c = y.eqn.res';
    case(1)
        c(1) = x(1)*x(2)*x(3)*x(4);
        c(2) = x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2;
    case(2)
        c(1) = sum(x);
    case(3)
        c(1) = 0.1*x(1) - x(2);
    case(4)
        c(1) = x(1)*x(2);
        c(2) = x(1)^2 + x(2)^2;
    case(5)
        c(1) = 0.1*x(1)-x(2);
        c(2) = -10*x(1)+x(2);
    case(6)
        c(1) = 2 * x(1) + x(2);
        c(2) = x(1) + 2* x(2);
    case(7)
        c(1) = x(1) + x(2);
    case(8)
        c(1) = 9 - x(1)^2 - x(2);
    case(9)
        c(1) = x(1)^2;
    case(10)
        c(1) = 0.1*x(1)-x(2);
end
