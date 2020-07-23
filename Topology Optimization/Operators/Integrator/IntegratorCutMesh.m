classdef IntegratorCutMesh < Integrator
    
    properties (Access = private)        
        backgroundCutMesh
        cutMeshOfSubCellLocal
        cellContainingSubcell
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.backgroundCutMesh     = cParams.backgroundCutMesh; 
            obj.cutMeshOfSubCellLocal = cParams.cutMeshOfSubCellLocal;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
        end
        
        function rhs = integrate(obj,fNodal)
            feMesh = obj.backgroundCutMesh;
            xGauss = obj.computeGaussPoints();
            rhsCellsCut = obj.computeElementalRHS(fNodal,feMesh,xGauss);
            rhsCells    = obj.assembleSubcellsInCells(rhsCellsCut);
            rhs = obj.assembleIntegrand(rhsCells);
        end
        
    end
    
    methods (Access = private)
        
        function xGauss = computeGaussPoints(obj)
            q = obj.computeQuadrature();
            m = obj.cutMeshOfSubCellLocal;
            xGauss = m.computeXgauss(q.posgp);
        end
        
        function rhsCells = assembleSubcellsInCells(obj,rhsCut)           
            nnode = size(obj.globalConnec,2);
            nelem = size(obj.globalConnec,1);
            cellNum = obj.cellContainingSubcell;
            totalInt = zeros(nelem,nnode);
            for iNode = 1:nnode
                int = rhsCut(:,iNode);
                intGlobal = accumarray(cellNum,int,[nelem,1],@sum,0);
                totalInt(:,iNode) = totalInt(:,iNode) + intGlobal;
            end
            rhsCells = totalInt;
        end
        
    end
    
end

