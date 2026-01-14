classdef NonLocalDamageFunctional < handle
    
    properties (Access = private)
        mesh
        cw
        Gc
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
            F   = ((obj.Gc*obj.l0)/(2*obj.cw))*int.compute(phiGradSquaredFun);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            J = IntegrateRHS(@(v) ((obj.Gc*obj.l0)/(obj.cw))*DP(Grad(v),Grad(phi)),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end
        
        function H = computeHessian(obj,~,quadOrder)
            H = IntegrateLHS(@(u,v) ((obj.Gc*obj.l0)/(obj.cw))*DP(Grad(v),Grad(u)),obj.testPhi,obj.testPhi,obj.mesh,'Domain',quadOrder);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;            
            obj.testPhi  = copy(cParams.testSpace.phi);
            obj.cw       = cParams.dissipation.constant;
            obj.Gc       = cParams.Gc;
            obj.l0       = cParams.l0;
        end
        
    end
end
