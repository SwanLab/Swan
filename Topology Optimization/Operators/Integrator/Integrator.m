classdef Integrator < handle

    properties (GetAccess = protected, SetAccess = protected)
       npnod 
       globalConnec
       mesh       
    end
    
    methods (Static, Access = public)
        
        function obj = create(cParams)
            obj = IntegratorFactory.create(cParams);
        end
        
    end

    methods (Access = protected)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.npnod              = cParams.npnod;
        end
        
        function quadrature = computeQuadrature(obj)
            quadrature = Quadrature.set(obj.mesh.type);
            quadrature.computeQuadrature('LINEAR');
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
        
        function rhsC = computeElementalRHS(obj,fNodal,xGauss,connec,type)
            s.fNodal         = fNodal;
            s.xGauss         = xGauss;
            s.mesh           = obj.mesh;
            s.connec         = connec;
            s.type           = type;
            rhs = RHSintegrator(s);
            rhsC = rhs.integrate();
        end        
        
        
    end
    
    
    
end

