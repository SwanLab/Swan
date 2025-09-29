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
        
        function F = computeCost(obj,phi,quadOrder)        
            phiGradSquaredFun = norm(Grad(phi)).^2;
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = (obj.constant*(obj.l0/2))*int.compute(phiGradSquaredFun);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            J = IntegrateRHS(@(v) (obj.constant*obj.l0)*DP(Grad(v),Grad(phi)),obj.testPhi,obj.mesh,quadOrder);
        end
        
        function H = computeHessian(obj,~,quadOrder)
            H = IntegrateLHS(@(u,v) (obj.constant*obj.l0)*DP(Grad(v),Grad(u)),obj.testPhi,obj.testPhi,obj.mesh,'Domain',quadOrder);
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