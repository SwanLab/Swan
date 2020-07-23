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
            connec   = obj.mesh.connec;     
            type     = obj.mesh.type;     
            xGauss   = obj.computeGaussPoints();
            rhsCells = obj.computeElementalRHS(fNodal,xGauss,connec,type);
            rhs = obj.assembleIntegrand(rhsCells);
        end
        
    end
    
    methods (Access = private)
        
        function xGauss = computeGaussPoints(obj)
            q = obj.computeQuadrature();
            xGauss = repmat(q.posgp,[1,1,obj.mesh.nelem]);
        end
        
    end
    
end