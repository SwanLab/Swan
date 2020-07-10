classdef IntegratorCutMesh < Integrator
    
    properties (Access = private)
        RHScells
        RHScellsCut
        backgroundMesh
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.backgroundMesh = cParams.meshBackground;
        end
        
        function rhs = integrate(obj,Fnodal)
            rhsCellsCut = obj.computeElementalCutRHS(Fnodal);
            rhsCells    = obj.assembleSubcellsInCells(rhsCellsCut);            
            rhs = obj.assembleIntegrand(rhsCells);
        end
        
    end
    
    methods (Access = private)
        
        function rhsV = computeElementalCutRHS(obj,fNodal)
            s.fNodal         = fNodal;
            s.xGauss         = obj.computeUnfittedGaussPoints();            
            s.quadrature     = obj.computeQuadrature(obj.mesh.geometryType);
            s.backgroundMesh = obj.backgroundMesh;
            s.mesh           = obj.mesh;            
            s.feMesh         = obj.mesh.cutMeshOfSubCellGlobal();
            rhs = RHSintegrator(s);
            rhsV = rhs.integrate();                  
        end
        
        function xGauss = computeUnfittedGaussPoints(obj)
            q = obj.computeQuadrature(obj.mesh.geometryType);
            m = obj.mesh.cutMeshOfSubCellLocal;
            xGauss = m.computeXgauss(q.posgp);
        end
        
        function rhsCells = assembleSubcellsInCells(obj,rhsCut)
            nnode = obj.backgroundMesh.nnode;
            nelem = obj.backgroundMesh.nelem;
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
            npnod  = obj.backgroundMesh.npnod;
            nnode  = obj.backgroundMesh.nnode;
            connec = obj.backgroundMesh.connec;
            f = zeros(npnod,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[npnod,1],@sum,0);
            end
        end
        
    end

end

