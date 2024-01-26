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

        function computeFirstStepLength(obj,g,Vtar,DM)
            phi0 =  obj.phi.fun.fValues;
            V0 = (g.value+1)*Vtar;
            tUp   = 1;
            tDown = 0;

            obj.tau = tUp;
            phiUp = obj.update(DM);
            obj.tau = tDown;
            phiDown = obj.update(DM);
            obj.phi.update(phiUp);
            g.computeFunctionAndGradient();
            Vup = (g.value+1)*Vtar;
            obj.phi.update(phiDown);
            g.computeFunctionAndGradient();
            VDown = (g.value+1)*Vtar;
            obj.phi.update(phi0);
            dVup = V0-Vup-0.05;
            dVdown = V0-VDown-0.05;

            hasFound = false;
            while not(hasFound)
                if obj.isFeasible(dVup)
                    tau0 = tUp;
                    hasFound = true;
                elseif obj.isFeasible(dVdown)
                    tau0 = tDown;
                    hasFound = true;
                else
                    tNew = (tUp+tDown)/2;
                    obj.tau = tNew;
                    phiNew = obj.update(DM);
                    obj.phi.update(phiNew);
                    g.computeFunctionAndGradient();
                    Vnew = (g.value+1)*Vtar;
                    DVnew = V0-Vnew-0.05;
                    if DVnew >=0
                        tUp = tNew;
                        dVup = DVnew;
                    else
                        tDown = tNew;
                        dVdown = DVnew;
                    end
                    obj.phi.update(phi0);
                end
            end
            g.computeFunctionAndGradient();
            obj.tau = tau0;
        end

            function is = isFeasible(obj,dV)
                cond1 = dV <= 0;
                cond2 = dV > -0.05;
                is = cond1 & cond2;
            end

            function is = isTooSmall(obj)
                is = obj.tau < 1e-10;
            end

        function increaseStepLength(obj,f)
            obj.tau = min(f*obj.tau,1);
        end

        function decreaseStepLength(obj)
            obj.tau = obj.tau/1.1;
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