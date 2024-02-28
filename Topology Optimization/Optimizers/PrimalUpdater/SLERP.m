classdef SLERP < handle

    properties (Access = public)
        tau
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        function obj = SLERP(cParams)
            obj.init(cParams);
        end

        function phi = update(obj,g,phi)     
            phiF   = obj.createP1Function(phi);
            gF     = obj.createP1Function(g);
            gN     = gF.normalize('L2');
            phiN   = phiF.normalize('L2');
            theta  = obj.computeTheta(phiN,gN);
            phiNew = obj.computeNewLevelSet(phiN,gN,theta);
            phi    = phiNew.fValues;
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
            obj.tau = obj.tau/1.1;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function f = createP1Function(obj,fV)
            s.mesh    = obj.mesh;
            s.fValues = fV;
            s.order   = 'P1';
            f         = LagrangianFunction(s);
        end

        function t = computeTheta(obj,phi,g)
            m = obj.mesh;
            phiG = ScalarProduct.computeL2(m,phi,g);
            t = max(acos(phiG),1e-14);
        end

        function pF = computeNewLevelSet(obj,phi,g,theta)
            k  = obj.tau;
            t  = theta;
            pN = phi.fValues;
            gN = g.fValues;
            a  = sin((1-k)*t)/sin(t);
            b  = sin(k*t)/sin(t);
            p  = a*pN + b*gN;
            pF  = obj.createP1Function(p);
        end

    end

end