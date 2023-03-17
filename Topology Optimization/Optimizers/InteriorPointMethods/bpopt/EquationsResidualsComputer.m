%% BPOPT Solver: Obtain equation residuals
classdef EquationsResidualsComputer < handle
    properties (Access = public)
        cRes
    end
    properties (Access = private)
        bp
        xC
    end

    methods (Access = public)
        function obj = EquationsResidualsComputer(cParams)
            obj.init(cParams);
        end
        function create(obj)
            obj.createCases()          
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
            obj.xC = cParams.x;
        end
        function createCases(obj)
            x = obj.xC;
            switch obj.bp.prob
                case(0)
                    %y = apm_details(bp.server,bp.app,x);
                    %obj.cRes = y.eqn.res';
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
end