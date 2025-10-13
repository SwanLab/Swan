classdef NeohookeanFunctional < handle

    properties (Access = private)
        lambda
        mu
        mesh
    end

    methods (Access = public)

        function obj = NeohookeanFunctional(cParams)
            obj.init(cParams)
        end

        function val = computeCost(obj,uFun,quadOrder)
            fun = obj.computeCostFunction(uFun);
            s.mesh = obj.mesh;
            s.quadType = quadOrder;
            int = IntegratorFunction(s);
            val = int.compute(fun);
        end

        function Fint = computeGradient(obj,uFun,quadOrder)
            PK1 = obj.computeFirstPiola(uFun);
            Fint = IntegrateRHS(@(v) DDP(Grad(v),PK1),uFun,obj.mesh,quadOrder);
        end

        function hess = computeHessian(obj,uFun,quadOrder)
            % This is the LINEALIZED hessian (Holzapfel, 401)
            % See  Holzapfel, 396
            Aneofun = obj.computeTangentConstitutive(uFun);
            hess = IntegrateLHS(@(u,v) DDP(Grad(v),DDP(Aneofun,Grad(u))),uFun,uFun,obj.mesh,quadOrder);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.lambda = cParams.material.lambda;
            obj.mu     = cParams.material.mu;
            obj.mesh   = cParams.mesh;
        end

        function fun = computeCostFunction(obj,uFun)
            nDim = obj.mesh.ndim;
            [~,F] = obj.computeDeformationGradient(uFun);
            C = F'*F;
            trC = trace(C);
            jac = Det(F);
            fun = obj.mu./2.*(trC - nDim) + ...
                  -obj.mu.*log(jac) + ...
                  obj.lambda./2.*(log(jac)).^2;
        end

        function PK1 = computeFirstPiola(obj,uFun)
            [~,F] = obj.computeDeformationGradient(uFun);
            invFt = Inv(F');
            jac = Det(F);
            PK1 = obj.mu.*(F-invFt) + obj.lambda.*log(jac).*invFt;
        end

        function Aneo = computeTangentConstitutive(obj,uFun)
            [I33,F] = obj.computeDeformationGradient(uFun);
            invF = Inv(F);
            jac = Det(F);
            Aneo = Expand(obj.lambda,4).*kronProd(invF',invF',[1 2 3 4]) + ...
                   Expand(obj.mu,4).*kronProd(I33,I33,[1 3 2 4]) + ...
                   Expand((obj.mu-obj.lambda.*log(jac)),4).*kronProd(invF',invF,[1 4 2 3]);
        end

        function [I33,F] = computeDeformationGradient(~,uFun)
            I33 = Identity(uFun);
            F = I33 + Grad(uFun)';
        end

    end

end