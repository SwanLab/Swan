classdef SLERP < handle

    properties (Access = public)
        tau
        Theta
        Alpha
        Beta
        boxConstraints
    end

    properties (Access = private)
        mesh
        volume
    end

    methods (Access = public)
        function obj = SLERP(cParams)
            obj.init(cParams);
        end

        function phi = update(obj,g,phi)  
            n         = obj.mesh.nnodes;
            phiF      = phi.levelSets;
            gF1       = obj.createP1Function(g(1:n));
            gF2       = obj.createP1Function(g(n+1:n*2));
            gF3       = obj.createP1Function(g(n*2+1:n*3));
            gN1       = gF1.normalize('L2');
            gN2       = gF2.normalize('L2');
            gN3       = gF3.normalize('L2');
            phiN1     = phiF{1,1}.fun.normalize('L2');
            phiN2     = phiF{1,2}.fun.normalize('L2');
            phiN3     = phiF{1,3}.fun.normalize('L2');
            gN        = obj.createP1Function([gN1.fValues,gN2.fValues,gN3.fValues]);
            phiN      = obj.createP1Function([phiN1.fValues,phiN2.fValues,phiN3.fValues]);
            theta1     = real(obj.computeTheta(phiN1,gN1));
            theta2     = real(obj.computeTheta(phiN2,gN2));
            theta3     = real(obj.computeTheta(phiN3,gN3));
            theta      = norm(norm(theta1,theta2),theta3);
            obj.Theta = theta;
            phiNew    = obj.computeNewLevelSet(phiN,gN,theta);
            phi.update(phiNew);
            obj.updateBoundsMultipliers(phi.fun);
        end

        function phi = updateFirstIter(obj,g,phi)
            phiF      = phi.fun;
            gF        = obj.createP1Function(g);
            gN        = gF.normalize('L2');
            phiN      = phiF.normalize('L2');
            theta     = obj.computeTheta(phiN,gN);
            phiNew    = obj.computeNewLevelSet(phiN,gN,theta);
            phi.update(phiNew);
        end

        function computeFirstStepLength(obj,g,ls,~)
            V0 = obj.volume.computeFunctionAndGradient(ls.levelSets{1,1});
            if abs(V0-1) <= 1e-10
                obj.computeLineSearchInBounds(g(1:obj.mesh.nnodes),ls.levelSets{1,1});
            else
                obj.tau = 1;
            end
        end

        function computeLineSearchInBounds(obj,g,ls)
            tLower  = 0;
            tUpper  = obj.computeInitialUpperLineSearch(g,ls);
            obj.tau = 0.5*(tUpper+tLower);
            V       = obj.computeVolumeFromTau(g,ls);
            delta   = abs(V-1);
            cond1   = delta<=1e-10;
            cond2   = delta>=0.05;
            while (cond1 || cond2)
                if cond1
                    tLower  = obj.tau;
                end
                if cond2
                    tUpper  = obj.tau;
                end
                obj.tau = 0.5*(tUpper+tLower);
                V       = obj.computeVolumeFromTau(g,ls);
                delta   = abs(V-1);
                cond1   = delta<=1e-10;
                cond2   = delta>=0.05;
            end
        end

        function tU = computeInitialUpperLineSearch(obj,g,ls)
            obj.tau = 1;
            V       = obj.computeVolumeFromTau(g,ls);
            delta   = abs(V-1);
            while delta<0.05
                obj.tau = obj.tau*2;
                V       = obj.computeVolumeFromTau(g,ls);
                delta   = abs(V-1);
            end
            tU = obj.tau;
        end

        function V = computeVolumeFromTau(obj,g,ls)
            lsAux = ls.copy();
            lsAux = obj.updateFirstIter(g,lsAux);
            V     = obj.volume.computeFunctionAndGradient(lsAux);
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
            obj.mesh = cParams.mesh;
            %obj.tau = cParams.tau;
            obj.createVolumeFunctional();
        end

        function createVolumeFunctional(obj)
            s.mesh         = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.volume     = VolumeFunctional(s);
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

        function p = computeNewLevelSet(obj,phi,g,theta)
            k  = obj.tau;
            t  = theta;
            pN = phi.fValues;
            gN = g.fValues;
            a  = sin((1-k)*t)/sin(t);
            b  = sin(k*t)/sin(t);
            p  = a*pN + b*gN;
            obj.Alpha = a;
            obj.Beta  = b;
        end

        function updateBoundsMultipliers(obj,xF)
            x                         = xF.fValues;
            obj.boxConstraints.lUB    = zeros(size(x));
            obj.boxConstraints.lLB    = zeros(size(x));
            obj.boxConstraints.refTau = 1;
        end

    end

end