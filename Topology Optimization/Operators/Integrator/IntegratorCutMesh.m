classdef IntegratorCutMesh < Integrator
    
    properties (Access = private)
        globalConnec
        npnod
        geometryType
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
        end
        
        function rhs = integrate(obj,Fnodal)
            rhsCellsCut = obj.computeElementalCutRHS(Fnodal);
            rhsCells    = obj.assembleSubcellsInCells(rhsCellsCut);
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
        
        function rhsV = computeElementalCutRHS(obj,fNodal)
            s.fNodal         = fNodal;
            s.xGauss         = obj.computeGaussPoints();
            s.quadrature     = obj.computeQuadrature(obj.mesh.geometryType);
            s.geometryType   = obj.mesh.cutMeshOfSubCellGlobal.geometryType;
            s.mesh           = obj.mesh;
            s.feMesh         = obj.mesh.cutMeshOfSubCellGlobal();
            rhs = RHSintegrator(s);
            rhsV = rhs.integrate();
        end
        
        function xGauss = computeGaussPoints(obj)
            q = obj.computeQuadrature(obj.mesh.geometryType);
            m = obj.mesh.cutMeshOfSubCellLocal;
            xGauss = m.computeXgauss(q.posgp);
        end
        
        function rhsCells = assembleSubcellsInCells(obj,rhsCut)           
            nnode = size(obj.globalConnec,2);
            nelem = size(obj.globalConnec,1);
            cellNum = obj.mesh.cellContainingSubcell;
            totalInt = zeros(nelem,nnode);
            for iNode = 1:nnode
                int = rhsCut(:,iNode);
                intGlobal = accumarray(cellNum,int,[nelem,1],@sum,0);
                totalInt(:,iNode) = totalInt(:,iNode) + intGlobal;
            end
            rhsCells = totalInt;
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

