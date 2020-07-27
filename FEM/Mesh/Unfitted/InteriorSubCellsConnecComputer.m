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
        nSubCellsByElem
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
            obj.computeNsubCellsByElem();
            obj.computeAllSubCellsConnec();
            obj.computeIsSubCellsInterior();
            obj.computeConnec();
            obj.computeXcoordsIso();
            obj.computeCellContainingSubCell();
        end
        
        function computeNsubCellsByElem(obj)
            switch size(obj.isSubCellInteriorParams.subCellCases,2)
                case 3
                    obj.nSubCellsByElem   = 3;            
                case 4
                    obj.nSubCellsByElem   = 4;
            end            
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
            s.nSubCellsByElem = obj.nSubCellsByElem; 
            subCellInt = IsSubCellInteriorComputer(s);
            subCellInt.compute();
            obj.isSubCellInterior = subCellInt.isSubCellInterior;
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
            localCell  = obj.computeLocalCellContainingSubCell();
            globalCell = localToGlobalCut(localCell); 
            obj.cellContainingSubcell = globalCell;
        end
        
        function cells = computeLocalCellContainingSubCell(obj)
            isInterior = obj.isSubCellInterior;            
            nCutElem = size(isInterior,2);
            localSubCellsInCell = repmat(1:nCutElem,obj.nSubCellsByElem,1);
            cells = localSubCellsInCell(isInterior);            
        end
        
                
    end
    
end
