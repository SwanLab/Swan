%% BPOPT Solver: Objective function value
function [ob] = bp_obj(bp,x)

switch bp.prob
    case(0)
        y = apm_details(bp.server,bp.app,x);
        ob = y.obj;
    case(1)
        ob = x(1)*x(4)*(x(1) + x(2) + x(3))  +  x(3);
    case(2)
        ob = 0;
        n = size(x,2);
        for i = 1:n,
            ob = ob + (x(i)-i)^2;
        end
    case(3)
        ob = x(1)^2 - 2*x(1)*x(2) + 4*x(2)^2;
    case(4)
        ob = x(2)*(5+x(1));
    case(5)
        ob = x(1)^2-2*x(1)*x(2)+4*x(2)^2;
    case(6)
        ob = x(1)^2 + 2 * x(2)^2;
    case(7)
        ob = (x(1)-5)^2 + (x(2)-5)^2;
    case(8)
        ob = (x(1)-5)^2;
    case(9)
        ob = (x(1)-5)^2;
    case(10)
        ob = x(1)^2-2*x(1)*x(2)+4*x(2)^2;
end