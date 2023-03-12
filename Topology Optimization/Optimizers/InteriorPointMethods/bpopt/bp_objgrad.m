%% BPOPT Solver: Obtain objective function gradient
classdef bp_objgrad < handle
    properties (Access = public)
        objGradient
    end
    properties (Access = private)
        xC 
        s
    end 

    methods (Access = public)
        function obj = bp_objgrad(cParams)
            obj.init(cParams);
        end
        function create(obj)
            obj.createObjectiveGradient();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.xC = cParams.x;
            obj.s = cParams.s;
        end
        function createObjectiveGradient(obj)
            x = obj.xC;
            switch bp.prob
            case(0)
                y = apm_details(bp.server,bp.app,x);
                g = y.obj_grad';        
            case(1)
                g(1) = x(4)*(2*x(1)+x(2)+x(3));
                g(2) = x(1)*x(4);
                g(3) = x(1)*x(4) + 1;
                g(4) = x(1)*(x(1)+x(2)+x(3));
            case(2)
                n = size(x,2);
                for i = 1:n
                    g(i) = 2*(x(i)-i);
                end
            case(3)
                g(1) = 2*x(1)-2*x(2);
                g(2) = -2*x(1)+8*x(2);
            case(4)
                g(1) = x(2);
                g(2) = (5+x(1));
            case(5)
                g(1) =  2*x(1)-2*x(2);
                g(2) = -2*x(1)+8*x(2);
            case(6)
                g(1) = 2.0 * x(1);
                g(2) = 4.0 * x(2);
            case(7)
                g(1) = 2 * (x(1)-5);
                g(2) = 2 * (x(2)-5);
            case(8)
                g(1) = 2 * (x(1) - 5);
                g(2) = 0;
            case(9)
                g(1) = 2 * (x(1) - 5);
            case(10)
                g(1) =  2*x(1)-2*x(2);
                g(2) = -2*x(1)+8*x(2);

        ns = size(obj.s,2);
        g = [g zeros(1,ns)];
        obj.objGradient = g;
        end
    end
end