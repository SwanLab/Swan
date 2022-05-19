classdef SLERP < handle

    properties (Access = public)
        tau
    end

    properties (Access = private)
        phi
        theta
    end

    methods (Access = public)

        function obj = SLERP(cParams)
            obj.init(cParams);
        end

        function x = update(obj,g,~)
            obj.computeTheta(g);
            x = obj.computeNewLevelSet(g);
        end

        function computeFirstStepLength(obj,~,~,~)
            obj.tau = 1;
        end

        function is = isTooSmall(obj)
            is = obj.tau < 1e-10;
        end

        function increaseStepLength(obj,f)
            obj.tau = min(f*obj.tau,1);
        end

        function decreaseStepLength(obj)
            obj.tau = obj.tau/2;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.phi = cParams.designVar;
        end

        function computeTheta(obj,g)
            p         = obj.phi.value;
            obj.theta = acos((g'*p)/(norm(g)*norm(p)));
        end

        function p = computeNewLevelSet(obj,g)
            k = obj.tau;
            t = obj.theta;
            p = obj.phi.value;
            a = sin((1-k)*t)*p;
            b = sin(k*t)*g/norm(g);
            p = (a + b)/sin(t); 
        end

    end

end