classdef NonLocalDamageFunctional < handle
    
    properties (Access = private)
        mesh
        %constant
        l0
        testPhi
    end
    
    methods (Access = public)
        
        function obj = NonLocalDamageFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunctional(obj,phi,quadOrder)        
            phiGradSquaredFun = norm(Grad(phi),2)^2;
            int = Integrator.create('Function',obj.mesh,quadOrder);
            %F   = (obj.constant*(obj.l0/2))*int.compute(phiGradSquaredFun);
            constant = updateConstant(obj,phi);
            F   = ((obj.l0/2))*int.compute(phiGradSquaredFun.*constant);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            RHS = RHSIntegrator.create(s);
            %J = (obj.constant*obj.l0)*RHS.compute(Grad(phi), obj.testPhi);
            constant = updateConstant(obj,phi);
            J = (obj.l0)*RHS.compute(Grad(phi).*constant, obj.testPhi);
        end
        
        function H = computeHessian(obj,phi,quadOrder)
            s.trial = obj.testPhi;
            s.test  = obj.testPhi;
            s.quadratureOrder = quadOrder;
            s.mesh = obj.mesh;
            %s.type = 'StiffnessMatrix';
            %LHS = LHSIntegrator.create(s);
            %H = (obj.constant*obj.l0)*LHS.compute();
            s.type = 'StiffnessMatrixWithFunction';
            s.function = obj.updateConstant(phi);
            LHS = LHSIntegrator.create(s);
            H = (obj.l0)*LHS.compute();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;            
            obj.testPhi  = copy(cParams.testSpace.phi);
            %obj.constant = cParams.dissipation.constant;
            obj.l0       = cParams.l0;
        end

        function c = updateConstant(obj,phi)
            isOverThreshold = (phi.fun>=0.75);
            c = (isOverThreshold).*(1.5*1.5*0.1) + (~isOverThreshold).*(1.5*0.1);
        end
        
    end
end