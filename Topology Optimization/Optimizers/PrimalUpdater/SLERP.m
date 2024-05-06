classdef SLERP < handle

    properties (Access = public)
        tau
        Theta
        boxConstraints
    end

    properties (Access = private)
        mesh
        volume
        filter
    end

    methods (Access = public)
        function obj = SLERP(cParams)
            obj.init(cParams);
            obj.createFilter();
        end

        function phi = update(obj,g,phi)   
            y = obj.computeRegularizedDensity(phi);
            phiF   = phi.fun;
            gF     = obj.createP1Function(g);
            gN     = gF.normalize('L2');
            phiN   = phiF.normalize('L2');
            theta  = obj.computeTheta(phiN,gN);
            obj.Theta = theta;
            phiNew = obj.computeNewLevelSet(phiN,gN,theta);
            phi.update(phiNew);
            x = obj.computeRegularizedDensity(phi);
            obj.updateBoundsMultipliers(x,y,g,phiNew);
        end

        function computeFirstStepLength(obj,g,ls,~)
            V0 = obj.volume.computeFunctionAndGradient(ls);
            if abs(V0-1) <= 1e-10
                obj.computeLineSearchInBounds(g,ls);
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
            lsAux = obj.update(g,lsAux);
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
            obj.createVolumeFunctional();
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.filter   = Filter.create(s);
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
        end

        function rhoe = computeRegularizedDensity(obj,phi)
            charFun = phi.obtainDomainFunction();
            rhoe    = obj.filter.compute(charFun,'QUADRATIC');
        end

        function updateBoundsMultipliers(obj,xF,yF,g,phi)
            x       = xF.fValues;
            y       = yF.fValues;
            t       = sum(abs(y-x))/sum(abs(g));
            isUBAct = phi<0 & g<0;
            isLBAct = phi>0 & g>0;
            lUB     = y-t*g-x;
            lLB     = x+t*g-y;

            lUB(~isUBAct | lUB<0)     = 0;
            lLB(~isLBAct | lLB<0)     = 0;
            obj.boxConstraints.lUB    = 0; % lUB
            obj.boxConstraints.lLB    = 0; % lLB
            obj.boxConstraints.refTau = t;
        end

    end

end