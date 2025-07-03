classdef NonLocalDamageFunctional < handle
    
    properties (Access = private)
        mesh
        constant
        l0
        testPhi
    end
    
    methods (Access = public)
        
        function obj = NonLocalDamageFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunctional(obj,phi,quadOrder)        
            phiGradSquaredFun = norm(Expand(Grad(phi)),2).^2;
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = (obj.constant*(obj.l0/2))*int.compute(phiGradSquaredFun);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            RHS = RHSIntegrator.create(s);
            J = (obj.constant*obj.l0)*RHS.compute(Grad(phi), obj.testPhi);
        end
        
        function H = computeHessian(obj,~,quadOrder)
            s.trial = obj.testPhi;
            s.test  = obj.testPhi;
            s.quadratureOrder = quadOrder;
            s.mesh = obj.mesh;
            s.type = 'StiffnessMatrix';
            LHS = LHSIntegrator.create(s);
            H = (obj.constant*obj.l0)*LHS.compute();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;            
            obj.testPhi  = copy(cParams.testSpace.phi);
            obj.constant = cParams.dissipation.constant;
            obj.l0       = cParams.l0;
        end
        
    end
end