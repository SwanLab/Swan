classdef Integrator_Simple < Integrator
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            obj.globalConnec = cParams.globalConnec;            
        end
        
        function LHS = computeLHS(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.globalConnec;
            s.npnod        = obj.npnod;
            lhs = LHSintegrator(s);
            LHS = lhs.compute();
        end
        
        function rhs = integrate(obj,fNodal)
            quadOrder = 'LINEAR';
            rhs = obj.integrateFnodal(fNodal,quadOrder);
        end
        
        function rhs = integrateFnodal(obj,fNodal,quadOrder)
            connec   = obj.mesh.connec;               
            xGauss   = obj.computeGaussPoints(quadOrder);
            type     = obj.mesh.type;                 
            fGauss   = obj.computeFgauss(fNodal,xGauss,connec,type);            
            rhs = obj.integrateFgauss(fGauss,xGauss,quadOrder);            
        end
        
        function rhs = integrateFgauss(obj,fGauss,xGauss,quadOrder)
            type     = obj.mesh.type;                 
            rhsCells = obj.computeElementalRHS(fGauss,xGauss,type,quadOrder);
            rhs = obj.assembleIntegrand(rhsCells);                        
        end
        
    end
    
    methods (Access = private)
        
        function xGauss = computeGaussPoints(obj,quadOrder)
            q = obj.computeQuadrature(quadOrder);
            xGauss = repmat(q.posgp,[1,1,obj.mesh.nelem]);
        end
        
    end
    
end