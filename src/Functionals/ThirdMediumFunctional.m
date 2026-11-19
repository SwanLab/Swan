classdef ThirdMediumFunctional < handle

    properties (Access = private)
        mesh
        density
        solidMu
        solidLambda
        gamma
        initialGamma
        finalGamma
        gammaReduction
        regularization
        neohookean
        regularizationHessian
    end

    methods (Access = public)

        function obj = ThirdMediumFunctional(cParams)
            obj.init(cParams);
            obj.gamma = obj.initialGamma;
            obj.createNeohookean();
        end

        function val = computeCost(obj,uFun,quadOrder)
            val = obj.neohookean.computeCost(uFun,quadOrder) + obj.computeRegularizationCost(uFun);
        end

        function Fint = computeGradient(obj,uFun,quadOrder)
            Fint = obj.neohookean.computeGradient(uFun,quadOrder) + obj.computeRegularizationGradient(uFun);
        end

        function hess = computeHessian(obj,uFun,quadOrder)
            K    = obj.computeRegularizationHessian(uFun);
            hess = obj.neohookean.computeHessian(uFun,quadOrder) + K;
        end

        function resetParameters(~)
        end

        function itHas = updateParameters(obj,~)
            itHas = obj.gamma <= obj.finalGamma;
            if ~itHas
                obj.gamma = max(obj.gamma/obj.gammaReduction,obj.finalGamma);
                obj.createNeohookean();
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            m = cParams.material;
            obj.mesh           = cParams.mesh;
            obj.density        = m.density;
            obj.solidMu        = m.mu;
            obj.solidLambda    = m.lambda;
            obj.initialGamma   = m.initialGamma;
            obj.finalGamma     = m.finalGamma;
            obj.gammaReduction = m.gammaReduction;
            obj.regularization = m.regularization;
        end

        function createNeohookean(obj)
            N      = obj.mesh.ndim;
            rho    = obj.density;
            muB    = obj.solidMu;
            kB     = obj.solidLambda + 2*muB/N;
            muA    = obj.gamma*muB;
            kA     = obj.gamma*kB;
            mu     = SimpAllInterpolator.computeMu(muA,muB,kA,kB,rho,N);
            kappa  = SimpAllInterpolator.computeKappa(muA,muB,kA,kB,rho,N);
            s.material.mu     = mu;
            s.material.lambda = kappa - (2/N)*mu;
            s.mesh            = obj.mesh;
            obj.neohookean    = NeohookeanFunctional(s);
        end

        function c = computeRegularizationCost(obj,uFun)
            v = obj.toDofVector(uFun);
            c = 0.5*v'*obj.computeRegularizationHessian(uFun)*v;
        end

        function g = computeRegularizationGradient(obj,uFun)
            g = obj.computeRegularizationHessian(uFun)*obj.toDofVector(uFun);
        end

        function K = computeRegularizationHessian(obj,uFun)
            if isempty(obj.regularizationHessian)
                w = obj.computeRegularizationWeight();
                f = @(u,v) w.*DDP(Hess(v),Hess(u));
                obj.regularizationHessian = IntegrateLHS(f,uFun,uFun,obj.mesh,'Domain',2);
            end
            K = obj.regularizationHessian;
        end

        function w = computeRegularizationWeight(obj)
            k = obj.regularization;
            w = k - k.*obj.density;
        end

    end

    methods (Static, Access = private)

        function v = toDofVector(uFun)
            v = reshape(uFun.fValues',[],1);
        end

    end

end
