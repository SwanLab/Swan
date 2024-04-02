classdef ProjectedGradient < handle

    properties (Access = public)
        tau
    end

    properties (Access = private)
        upperBound
        lowerBound
    end
    
    methods (Access = public)

        function obj = ProjectedGradient(cParams)
            obj.init(cParams);
        end

        function x = update(obj,g,x)
            ub = obj.upperBound;
            lb = obj.lowerBound;
            t  = obj.tau;
            x  = x - t*g;
            x  = min(ub,max(x,lb));
        end

        function computeFirstStepLength(obj,g,x,f)
            xVal    = x.fun.fValues;
            obj.tau = f*sqrt(norm(g)/norm(xVal));
        end
        
        function increaseStepLength(obj,f)
            obj.tau = f*obj.tau;
        end

        function decreaseStepLength(obj)
            obj.tau = obj.tau/2;
        end

        function is = isTooSmall(obj)
            is = obj.tau < 1e-10;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.upperBound     = cParams.ub;
            obj.lowerBound     = cParams.lb;
        end

    end

end