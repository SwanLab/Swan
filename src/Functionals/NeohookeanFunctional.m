classdef NeohookeanFunctional < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        lambda
        mu
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = NeohookeanFunctional(cParams)
            obj.init(cParams)
        end

        function val = computeCost(obj, uFun)      
            nDim = obj.mesh.ndim;
            [~,F] = obj.computeDeformationGradient(uFun);
            C = F'*F;
            trC = trace(C);
            jac = Det(F);
            s.quadType = 2;
            s.mesh = obj.mesh;           
            int = IntegratorFunction(s);
            f = obj.mu./2.*(trC - nDim) - obj.mu.*log(jac) + obj.lambda./2.*(log(jac)).^2;
            val = int.compute(f);
        end
        
        function Fint = computeGradient(obj, uFun)
            PK1 = obj.computeFirstPiola(uFun);
            s.mesh = obj.mesh;
            s.quadratureOrder = 3;
            s.type = 'ShapeDerivativeTensor';
            rhs = RHSIntegrator.create(s);
            Fint = rhs.compute(PK1,uFun); 
        end

        function hess = computeHessian(obj, uFun)
            % This is the LINEALIZED hessian (Holzapfel, 401)
            % See  Holzapfel, 396
            Aneofun = obj.computeTangentConstitutive(uFun);
            s.quadratureOrder = 3;
            s.test  = uFun;
            s.trial = uFun;
            s.mesh  = obj.mesh;
            s.type ='StiffnessFiniteStrain';
            lhs = LHSIntegrator.create(s);
            hess = lhs.compute(Aneofun);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.lambda = cParams.material.lambda;
            obj.mu     = cParams.material.mu;
            obj.mesh   = cParams.mesh;
        end

        function PK1 = computeFirstPiola(obj,uFun,xG)
            [~,F] = obj.computeDeformationGradient(uFun);
             invFt = Inv(F');
             jac = Det(F);
             PK1 = obj.mu.*(F-invFt) + obj.lambda.*log(jac).*invFt;
        end

        function Aneo = computeTangentConstitutive(obj,uFun)
            [I33,F] = obj.computeDeformationGradient(uFun);
            invF = Inv(F);
            jac = Det(F);
            Aneo = Expand(obj.lambda,4).*OP(invF', invF') + ...
                   Expand(obj.mu,4).*kronTop(I33,I33) + ...
                   Expand((obj.mu-obj.lambda.*log(jac)),4).*kronBot(invF', invF);
        end

        function [I33,F] = computeDeformationGradient(obj, uFun)
            I33 = Identity(uFun);
            F = I33 + Grad(uFun)';
        end     
    end
    
end