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
            phiN              = obj.normalizeLevelSets(phi);
            gN                = obj.createNormalizedGradient(phi,g);
            theta             = obj.computeThetaNorm(phi,phiN,gN);
            obj.Theta         = theta;
            [phiNvals,gNvals] = obj.computePhiAndGradientValues(phi,phiN,gN);
            phiNew            = obj.computeNewLevelSet(phiNvals,gNvals,theta);
            phi.update(phiNew);
            obj.updateBoundsMultipliers(phi.fun);
        end

        function computeFirstStepLength(obj,g,ls,~)
            switch class(ls)
                case 'LevelSet'
                    gClass  = g;
                    lsClass = ls;
                case 'MultiLevelSet'
                    gClass  = g(1:obj.mesh.nnodes);
                    lsClass = ls.levelSets{1,1};
            end
            V0 = obj.volume.computeFunctionAndGradient(lsClass);
            if abs(V0-1) <= 1e-10
                obj.computeLineSearchInBounds(gClass,lsClass);
            else
                obj.tau = 1;
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
            obj.createVolumeFunctional();
        end

        function createVolumeFunctional(obj)
            s.mesh         = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.volume     = VolumeFunctional(s);
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

        function phiN = normalizeLevelSets(obj,phi)
            switch class(phi)
                case 'LevelSet'
                    phiF = phi.fun;
                    phiN = phiF.normalize('L2');
                case 'MultiLevelSet'
                    phiF    = phi.levelSets;
                    phiN{1} = phiF{1,1}.fun.normalize('L2');
                    phiN{2} = phiF{1,2}.fun.normalize('L2');
                    phiN{3} = phiF{1,3}.fun.normalize('L2');
            end
        end

        function fN = createNormalizedGradient(obj,ls,fV)
            s.mesh  = obj.mesh;
            s.order = 'P1';
            switch class(ls)
                case 'LevelSet'
                    s.fValues = fV;
                    f         = LagrangianFunction(s);
                    fN        = f.normalize('L2');
                case 'MultiLevelSet'
                    nLS = length(ls.levelSets);
                    fV  = reshape(fV,[],nLS);
                    for i = 1:nLS
                        s.fValues = fV(:,i);
                        f         = LagrangianFunction(s);
                        fN{i}     = f.normalize('L2');
                    end
            end
        end

        function t = computeTheta(obj,phi,g)
            m    = obj.mesh;
            phiG = ScalarProduct.computeL2(m,phi,g);
            t    = max(acos(phiG),1e-14);
        end

        function t = computeThetaNorm(obj,phi,phiN,gN)
            switch class(phi)
                case 'LevelSet'
                    t = obj.computeTheta(phiN,gN);
                case 'MultiLevelSet'
                    nLS = length(phi.levelSets);
                    t   = obj.computeTheta(phiN{1},gN{1});
                    for i = 2:nLS
                        ti = obj.computeTheta(phiN{i},gN{i});
                        t  = norm(t,ti);
                    end
            end
        end

        function [phiNv,gNv] = computePhiAndGradientValues(obj,phi,phiN,gN)
            switch class(phi)
                case 'LevelSet'
                    phiNv = phiN.fValues;
                    gNv   = gN.fValues;
                case 'MultiLevelSet'
                    nLS   = length(phi.levelSets);
                    phiNv = [];
                    gNv   = [];
                    for i = 1:nLS
                        phiNv = [phiNv;phiN{i}.fValues];
                        gNv   = [gNv;gN{i}.fValues];
                    end
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

end