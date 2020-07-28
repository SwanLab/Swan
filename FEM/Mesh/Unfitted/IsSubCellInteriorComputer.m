classdef IsSubCellInteriorComputer < handle
    
    properties (Access = public)
        isSubCellInterior
        localCellContainingSubCell
    end
    
    properties (Access = private)
        subCellCases
        levelSet
        allNodesInElem
        nSubCellsByElem        
    end
    
    methods (Access = public)
        
        function obj = IsSubCellInteriorComputer(cParams)
            obj.init(cParams);
            obj.computeNsubCellsByElem();
        end
        
        function compute(obj)
            obj.computeIsSubCellInterior();
            obj.computeLocalCellContainingSubCell()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.subCellCases    = cParams.subCellCases;
            obj.levelSet        = cParams.levelSet;
            obj.allNodesInElem  = cParams.allNodesInElem;
        end
        
        function computeNsubCellsByElem(obj)
            switch size(obj.subCellCases,2)
                case 3
                    obj.nSubCellsByElem   = 3;            
                case 4
                    obj.nSubCellsByElem   = 4;
            end            
        end                
        
        function computeIsSubCellInterior(obj)
            nElem   = size(obj.subCellCases,1);
            isoNode = obj.computeIsoNode();
            isTriInt  = obj.computeIsSubCellTriangleInterior(isoNode);
            nSubCells = obj.nSubCellsByElem;
            itIs = false(nSubCells,nElem);
            itIs(1,isTriInt)  = true;
            itIs(2,~isTriInt) = true;
            itIs(3,~isTriInt) = true;
            obj.isSubCellInterior = itIs;            
        end
        
        function isoNode = computeIsoNode(obj)
            nCases  = size(obj.subCellCases,2);
            nElem   = size(obj.subCellCases,1);
            nodeIso = zeros(nElem,1);
            for icase = 1:nCases
                subCells = obj.subCellCases(:,icase);
                nodeIso(subCells,1) = icase;
            end
            allElems(:,1) = 1:nElem;
            t = sub2ind([nElem obj.nSubCellsByElem],allElems,nodeIso);
            isoNode = obj.allNodesInElem(t);
        end
        
        function itIs = computeIsSubCellTriangleInterior(obj,node)
            isoNodeIsFull = obj.levelSet(node) < 0;
            itIs = isoNodeIsFull;
        end
        
        function computeLocalCellContainingSubCell(obj)
            isInterior = obj.isSubCellInterior;    
            nCutElem   = size(isInterior,2);
            localSubCellsInCell = repmat(1:nCutElem,obj.nSubCellsByElem,1);
            cells = localSubCellsInCell(isInterior);         
            obj.localCellContainingSubCell = cells;
        end
        
    end
    
    
end