classdef ProjectedGradient < handle

    properties (Access = public)
        tau
        lUB
        lLB
    end

    properties (Access = private)
        upperBound
        lowerBound
    end
    
    methods (Access = public)
        function obj = ProjectedGradient(cParams)
            obj.init(cParams);
        end

        function x = update(obj,g,y)
            ub = obj.upperBound;
            lb = obj.lowerBound;
            t  = obj.tau;
            y  = y - t*g;
            x  = min(ub,max(y,lb));
            obj.updateBoundsMultipliers(x,y);
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
            obj.upperBound = cParams.ub;
            obj.lowerBound = cParams.lb;
        end

        function updateBoundsMultipliers(obj,x,y)
            dyx            = y-x;
            dxy            = x-y;
            obj.lUB        = zeros(size(x));
            obj.lLB        = zeros(size(x));
            obj.lUB(dyx>0) = dyx(dyx>0);
            obj.lLB(dxy>0) = dxy(dxy>0);
        end
    end

end