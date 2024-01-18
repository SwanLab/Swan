classdef SLERP < handle

    properties (Access = public)
        tau
    end

    properties (Access = private)
        phi
        theta
        epsilon
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
            obj.tau = 0.961;
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
            obj.phi     = cParams.designVar;
            obj.epsilon = cParams.uncOptimizerSettings.scalarProductSettings.femSettings.epsilon;
        end

        function computeTheta(obj,g)
            m         = obj.phi.mesh;            
            pN        = obj.normalizeFunction(obj.phi.fun.fValues);
            gN        = obj.normalizeFunction(g);
            s.fValues = pN;
            s.mesh    = obj.phi.mesh;
            pNfun     = P1Function(s);
            s.fValues = gN;
            gNfun     = P1Function(s);
            phiG      = ScalarProduct.computeH1(m,pNfun,gNfun,obj.epsilon);
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
            m         = obj.phi.mesh;
            s.fValues = x;
            s.mesh    = m;
            xFun      = P1Function(s);
            norm      = Norm.computeH1(m,xFun,obj.epsilon);
            xNorm     = sqrt(norm);
            x         = x/xNorm;
        end
    end

end