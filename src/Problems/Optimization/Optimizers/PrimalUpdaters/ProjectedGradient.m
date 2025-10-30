classdef ProjectedGradient < handle

    properties (Access = public)
        tau
        boxConstraints
    end

    properties (Access = private)
        upperBound
        lowerBound
        tauMax
    end
    
    methods (Access = public)
        function obj = ProjectedGradient(cParams)
            obj.init(cParams);
        end

        function rho = update(obj,g,rho)  
            y  = rho.fun.fValues;
            ub = obj.upperBound;
            lb = obj.lowerBound;
            t  = obj.tau;
            y  = y - t*g;
            x  = min(ub,max(y,lb));
            obj.updateBoundsMultipliers(x,y);
            rho.update(x);
        end

        function computeFirstStepLength(obj,g,x,f)
            if size(x,2) == 2
                x = x{1};
            end
            xVal    = x.fun.fValues;
            obj.tau = min(f*sqrt(norm(g)/norm(xVal)),obj.tauMax);
        end
        
        function increaseStepLength(obj,f)
            obj.tau = min(f*obj.tau,obj.tauMax);
        end

        function decreaseStepLength(obj)
            obj.tau = obj.tau/2;
        end

        function is = isTooSmall(obj)
            is = obj.tau < 1e-10;
        end

        function updateBounds(obj,ub,lb)
            if ~isnumeric(ub)
                obj.upperBound = ub.fValues;
            else
                obj.upperBound = ub;
            end
            if ~isnumeric(lb)
                obj.lowerBound = lb.fValues;
            else
                obj.lowerBound = lb;
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.upperBound = cParams.ub;
            obj.lowerBound = cParams.lb;
            obj.tauMax     = cParams.tauMax;
            obj.tau        = cParams.tau;
        end

        function updateBoundsMultipliers(obj,x,y)
            t          = obj.tau;
            dyx        = y-x;
            dxy        = x-y;
            lUB        = zeros(size(x));
            lLB        = zeros(size(x));
            lUB(dyx>0) = dyx(dyx>0);
            lLB(dxy>0) = dxy(dxy>0);
            obj.boxConstraints.lUB    = lUB;
            obj.boxConstraints.lLB    = lLB;
            obj.boxConstraints.refTau = t;
        end
    end

end