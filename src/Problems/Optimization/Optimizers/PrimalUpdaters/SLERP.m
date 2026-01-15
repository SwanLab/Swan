
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
    end

    methods (Access = public)
        function obj = SLERP(cParams)
            obj.init(cParams);
        end

        function phi = update(obj,g,phi)
            ls                = phi.obtainVariableInCell();
            phiN              = obj.normalizeLevelSets(ls);
            gN                = obj.createNormalizedGradient(ls,g);
            theta             = obj.computeThetaNorm(phiN,gN);
            obj.Theta         = theta;
            [phiNvals,gNvals] = obj.computePhiAndGradientValues(phiN,gN);
            phiNew            = obj.computeNewLevelSet(phiNvals,gNvals,theta);
            phi.update(phiNew);
            obj.updateBoundsMultipliers(phi.fun);
        end

        function computeFirstStepLength(obj,g,ls,~)
            [lsClass,gClass] = obj.getLevelSetAndGradientForVolume(ls,g);
            V0 = lsClass.computeVolume();
            if abs(V0-1) <= 1e-10
                obj.computeLineSearchInBounds(gClass,lsClass);
            else
                obj.tau = 0.1;
            end
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
            lsAux = copy(ls);
            lsAux = obj.update(g,lsAux);
            V     = lsAux.computeVolume();
        end

        function fN = createNormalizedGradient(obj,ls,fV)
            s.mesh  = obj.mesh;
            s.order = 'P1';
            nLS     = length(ls);
            fV      = reshape(fV,[],nLS);
            fN      = cell(nLS,1);
            for i = 1:nLS
                s.fValues = fV(:,i);
                f         = LagrangianFunction(s);
                fN{i}     = f.normalize('L2');
            end
        end

        function t = computeTheta(obj,phi,g)
            phiG = ScalarProduct(phi,g,'L2');
            t    = max(acos(phiG),1e-14);
        end

        function t = computeThetaNorm(obj,phiN,gN)
            t = 0;
            for i = 1:length(phiN)
                ti = obj.computeTheta(phiN{i},gN{i});
                t  = norm([t,ti]);
            end
        end

        function p = computeNewLevelSet(obj,pN,gN,theta)
            k  = obj.tau;
            t  = theta;
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

    methods (Static, Access = private)
        function [lsClass,gClass] = getLevelSetAndGradientForVolume(ls,g)
            lsC     = ls.obtainVariableInCell();
            lsClass = lsC{1};
            n       = length(lsClass.fun.fValues);
            gClass  = g(1:n);
        end

        function phiN = normalizeLevelSets(ls)
            nLS   = length(ls);
            phiN  = cell(nLS,1);
            for i = 1:nLS
                phiN{i} = ls{i}.fun.normalize('L2');
            end
        end

        function [phiNv,gNv] = computePhiAndGradientValues(phiN,gN)
            phiNv = [];
            gNv   = [];
            for i = 1:length(phiN)
                phiNv = [phiNv;phiN{i}.fValues];
                gNv   = [gNv;gN{i}.fValues];
            end
        end
    end

end