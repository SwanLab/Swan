classdef Integrator_Simple < Integrator
    
    properties (Access = private)
        globalConnec
        npnod

        
        RHScells
        RHSsubcells
        backgroundMesh        
    end
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.npnod          = cParams.npnod;            
            obj.globalConnec   = cParams.globalConnec;
        end
        
        function LHS = computeLHS(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.globalConnec;
            s.npnod        = obj.npnod;
            lhs = LHSintegrator(s);
            LHS = lhs.compute();
        end
        
        function norm = computeL2norm(obj,f)
            fv(:,1) = f(:);
            M = obj.LHS;
            norm = fv'*M*fv;
        end
        
        function RHS = integrate(obj,F)
            obj.computeElementalRHS(F);
            RHS = obj.assembleIntegrand();
        end
        
    end
    
    methods (Access = private)
        
        function computeElementalRHS(obj,fNodal)
            s.fNodal         = fNodal;
            s.xGauss         = obj.computeGaussPoints();   
            s.quadrature     = obj.computeQuadrature(obj.mesh.geometryType);            
            s.geometryType   = obj.backgroundMesh.geometryType;
            s.mesh           = obj.mesh;
            s.feMesh         = obj.mesh;
            rhs = RHSintegrator(s);
            obj.RHScells = rhs.integrate();        
        end
        
        function xGauss = computeGaussPoints(obj)
            q = obj.computeQuadrature(obj.mesh.geometryType);
            m = obj.mesh;
            xGauss = repmat(q.posgp,[1,1,m.nelem]);
        end             
        
        function f = assembleIntegrand(obj)
            integrand = obj.RHScells;
            ndofs  = obj.backgroundMesh.npnod;
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