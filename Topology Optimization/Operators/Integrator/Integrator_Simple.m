classdef Integrator_Simple < Integrator
    
    properties (Access = private)
        globalConnec        
        npnod
        geometryType
    end
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
        end
        
        function LHS = computeLHS(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.globalConnec;
            s.npnod        = obj.npnod;
            lhs = LHSintegrator(s);
            LHS = lhs.compute();
        end
        
        function rhs = integrate(obj,F)
            rhsCells = obj.computeElementalRHS(F);
            rhs = obj.assembleIntegrand(rhsCells);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.globalConnec = cParams.globalConnec;            
            obj.npnod        = cParams.npnod;
            obj.geometryType = cParams.geometryType;
        end
        
        function rhsC = computeElementalRHS(obj,fNodal)
            s.fNodal         = fNodal;
            s.xGauss         = obj.computeGaussPoints();
            s.quadrature     = obj.computeQuadrature(obj.mesh.geometryType);
            s.geometryType   = obj.geometryType;
            s.mesh           = obj.mesh;
            s.feMesh         = obj.mesh;
            rhs = RHSintegrator(s);
            rhsC = rhs.integrate();
        end
        
        function xGauss = computeGaussPoints(obj)
            m = obj.mesh;
            q = obj.computeQuadrature(m.geometryType);
            xGauss = repmat(q.posgp,[1,1,m.nelem]);
        end
        
        function f = assembleIntegrand(obj,rhsCells)
            integrand = rhsCells;
            ndofs  = obj.npnod;
            connec = obj.globalConnec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
        
    end
    
end