classdef IsSubCellInteriorComputer < handle
    
    properties (Access = public)
        isSubCellInterior
    end
    
    properties (Access = private)
        subCellCases
        levelSet
        allNodesInElem
        nSubCellsByElem
        nodesInSubCells
    end
    
    methods (Access = public)
        
        function obj = IsSubCellInteriorComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            nElem   = size(obj.subCellCases,1);
            isoNode = obj.computeIsoNode();
            isTriInt  = obj.computeIsSubCellTriangleInterior(isoNode);
            nSubCells = obj.nSubCellsByElem;
            itIs = false(nSubCells,nElem);
            itIs(1,isTriInt)  = true;
            itIs(2,~isTriInt) = true;
            itIs(3,~isTriInt) = true;
            obj.isSubCellInterior = itIs;
%             subCellInt = false(nSubCells,nElem);
%             for iSubCell = 1:3
%                 nodes = squeeze(obj.nodesInSubCells(:,iSubCell,:));
%                 ls = obj.levelSet(nodes);
%                 isT = any(ls);
%                 subCellInt(iSubCell,:) = isT;
%             end
%             
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.subCellCases    = cParams.subCellCases;
            obj.levelSet        = cParams.levelSet;
            obj.allNodesInElem  = cParams.allNodesInElem;
            obj.nSubCellsByElem = cParams.nSubCellsByElem;
            obj.nodesInSubCells = cParams.nodesInSubCells;
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
        
    end
    
end