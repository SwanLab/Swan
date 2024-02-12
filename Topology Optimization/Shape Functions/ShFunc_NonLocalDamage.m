classdef ShFunc_NonLocalDamage < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        constant
        l0
    end
    
    methods (Access = public)
        
        function obj = ShFunc_NonLocalDamage(cParams)
            obj.init(cParams)
            
        end
        function F = computeFunction(obj,phi,quad)        
            phiGrad = phi.computeGradient(quad);
            phiGradSquared = sum(phiGrad.fValues.^2);
            
            s.fValues = phiGradSquared;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            phiGradSquaredFun = FGaussDiscontinuousFunction(s);
            
            q.mesh = obj.mesh;
            q.quadType = quad.order;
            q.type = 'Function';
            int = Integrator.create(q);
            F = (obj.constant*(obj.l0/2))*int.compute(phiGradSquaredFun);
        end
        
        function J = computeGradient(obj,phi,quad)
            PhiGrad = phi.computeGradient(quad);
            test = P1Function.create(obj.mesh,1);
            
            s.quadratureOrder = quad.order;
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            RHS = RHSintegrator.create(s);
            J = (obj.constant*(obj.l0/2))*RHS.compute(PhiGrad, test);
        end
        
        function H = computeHessian(obj,quadOrder)
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.quadratureOrder = quadOrder;
            s.mesh = obj.mesh;
            s.type = 'StiffnessMatrix';
            LHS = LHSintegrator.create(s);
            H = (obj.constant*(obj.l0/2))*LHS.compute();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.constant = cParams.constant;
            obj.l0 = cParams.l0;
        end
        
    end
end