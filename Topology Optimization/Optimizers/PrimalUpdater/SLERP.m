classdef SLERP < handle

    properties (Access = public)
        tau
    end

    properties (Access = private)
        phi
        theta
        scalar_product
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
            obj.tau = obj.tau/1.5;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.phi            = cParams.designVar;
            obj.scalar_product = cParams.uncOptimizerSettings.scalarProductSettings;
        end

        function computeTheta(obj,g)
            pN   = obj.normalizeFunction(obj.phi.fun.fValues);
            gN   = obj.normalizeFunction(g);
            phiG = obj.scalar_product.computeSP(pN,gN);
            obj.theta = max(acos(phiG),1e-14);
        end

        function p = computeNewLevelSet(obj,g)
            k  = obj.tau;
            t  = obj.theta;
            pN = obj.normalizeFunction(obj.phi.fun.fValues);
            gN = obj.normalizeFunction(g);
            a  = sin((1-k)*t)*pN;
            b  = sin(k*t)*gN;
            p  = (a + b)/sin(t);
            p  = obj.normalizeFunction(p);
        end

        function x = normalizeFunction(obj,x)
            norm2 = obj.scalar_product.computeSP(x,x);
            xNorm = sqrt(norm2);
            x = x/xNorm;
        end

    end

end