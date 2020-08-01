classdef IntegratorCutMesh < Integrator
    
    properties (Access = private)        
        backgroundMeshType
        xCoordsIso
        cellContainingSubcell
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.globalConnec          = cParams.globalConnec;            
            obj.xCoordsIso            = cParams.xCoordsIso;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
            obj.backgroundMeshType    = cParams.backgroundMeshType;                        
        end
        
        function rhs = integrate(obj,fNodal)
            c = obj.computeSubCellConnec();
            t = obj.backgroundMeshType;
            xGauss = obj.computeGaussPoints();
            rhsCellsCut = obj.computeElementalRHS(fNodal,xGauss,c,t);
            rhsCells    = obj.assembleSubcellsInCells(rhsCellsCut);
            rhs = obj.assembleIntegrand(rhsCells);
        end
        
    end
    
    methods (Access = private)
        
        function xGauss = computeGaussPoints(obj)
            q = obj.computeQuadrature();            
            s.connec = obj.computeSubCellsLocalConnec();
            s.fNodes = obj.computeSubCellsLocalCoord();
            s.type   = obj.mesh.type;
            x = FeFunction(s);
            xGauss = x.interpolateFunction(q.posgp);
        end
        
        function c = computeSubCellsLocalCoord(obj)
            coord = obj.xCoordsIso; 
            nDim  = size(coord,1);            
            c = reshape(coord,nDim,[])';             
        end
        
        function lConnec = computeSubCellsLocalConnec(obj)
            coord = obj.xCoordsIso;
            nElem = size(coord,3);
            nNode = size(coord,2);
            lConnec = reshape(1:nElem*nNode,nNode,nElem)';            
        end
        
        function c = computeSubCellConnec(obj)
            cells = obj.cellContainingSubcell;
            c = obj.globalConnec(cells,:);        
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

