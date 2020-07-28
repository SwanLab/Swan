classdef InteriorSubCellsConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        connec
        xCoordsIso
        cellContainingSubcell
    end
    
    properties (Access = private)
        isSubCellInterior        
        allSubCellsConnec        
        xNodesInSubCells    
    end
    
    properties (Access = private)
        localCells
        allSubCellsConnecParams
        isSubCellInteriorParams
        cutElems
    end
    
    methods (Access = public)
        
        function obj = InteriorSubCellsConnecComputer(cParams)
            obj.init(cParams);
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)        
            obj.allSubCellsConnecParams = cParams.allSubCellsConnecParams;
            obj.isSubCellInteriorParams = cParams.isSubCellInteriorParams;
            obj.cutElems                = cParams.cutElems;
        end
        
        function compute(obj)
            obj.computeAllSubCellsConnec();
            obj.computeIsSubCellsInterior();
            obj.computeConnec();
            obj.computeXcoordsIso();
            obj.computeCellContainingSubCell();
        end
        
        function computeAllSubCellsConnec(obj)
            s = obj.allSubCellsConnecParams;
            a = AllSubCellsConnecComputer(s);
            a.compute();            
            obj.allSubCellsConnec = a.allSubCellsConnec;  
            obj.xNodesInSubCells  = a.xNodesInSubCells;
        end
           
        function computeIsSubCellsInterior(obj)
            s = obj.isSubCellInteriorParams;
            c = IsSubCellInteriorComputer(s);
            c.compute();
            obj.isSubCellInterior = c.isSubCellInterior;
            obj.localCells        = c.localCellContainingSubCell;
        end

        function computeConnec(obj)
            allConnec  = obj.allSubCellsConnec;
            isInterior = obj.isSubCellInterior(:);
            obj.connec = allConnec(isInterior,:);           
        end
                
        function computeXcoordsIso(obj)
            isInterior = obj.isSubCellInterior(:);            
            x = obj.xNodesInSubCells;
            nDim = size(x,3);
            nElem = sum(isInterior);
            nNode = size(x,1);
            obj.xCoordsIso = zeros(nDim,nNode,nElem);
            for idim = 1:nDim
                obj.xCoordsIso(idim,:,:) = x(:,isInterior,idim);
            end
        end
        
        function computeCellContainingSubCell(obj)
            localToGlobalCut = obj.cutElems;
            lCells  = obj.localCells;
            globalCell = localToGlobalCut(lCells); 
            obj.cellContainingSubcell = globalCell;
        end
        
  
    end
    
end
