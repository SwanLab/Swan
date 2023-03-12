%% BPOPT Solver: Obtain equation residuals
classdef bp_res_stub < handle
    properties (Access = public)
        cRes
    end
    properties (Access = private)
    end

    methods (Access = public)
        function obj = bp_res_stub(cParams)
            obj.init(cParams);
        end
        function create(obj)
            obj.createCases()          
        end
    end
    methods (Access = private)
        function init(obj,cParams)

        end
        function createCases(obj)
            switch cParams.prob
            case(0)
                y = apm_details(bp.server,bp.app,x);
                obj.cRes = y.eqn.res';
            case(1)
                obj.cRes(1) = x(1)*x(2)*x(3)*x(4);
                obj.cRes(2) = x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2;
            case(2)
                obj.cRes(1) = sum(x);
            case(3)
                obj.cRes(1) = 0.1*x(1) - x(2);
            case(4)
                obj.cRes(1) = x(1)*x(2);
                obj.cRes(2) = x(1)^2 + x(2)^2;
            case(5)
                obj.cRes(1) = 0.1*x(1)-x(2);
                obj.cRes(2) = -10*x(1)+x(2);
            case(6)
                obj.cRes(1) = 2 * x(1) + x(2);
                obj.cRes(2) = x(1) + 2* x(2);
            case(7)
                obj.cRes(1) = x(1) + x(2);
            case(8)
                obj.cRes(1) = 9 - x(1)^2 - x(2);
            case(9)
                obj.cRes(1) = x(1)^2;
            case(10)
                obj.cRes(1) = 0.1*x(1)-x(2);
        end
    end
end