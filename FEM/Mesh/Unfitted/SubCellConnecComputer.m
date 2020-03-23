classdef SubCellConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        connec
    end
    
    properties (Access = private)
        coord
        nElem
        connecTotal
        
        allSubCellsConnec
        
        elemCases
        levelSet
        isoNode
        vertexNodesInElem
        cutNodesInElem
        allNodesInElem
        
        
        nSubCellNodes
        nSubCells
        nSubCellsByElem
        nCases
        
        nSubCases
        nSubCellsByQuad
        
        nnodeT
        
        isSubCellInterior
    end
    
    methods (Access = public)
        
        function obj = SubCellConnecComputer(cParams)
            obj.init(cParams);
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.cutNodesInElem   = cParams.cutNodePerElemen;
            obj.vertexNodesInElem = cParams.backgroundConnec;
            obj.nElem            = cParams.nElem;
            obj.elemCases        = cParams.elemCases;
            obj.coord            = cParams.coord;
            obj.levelSet         = cParams.levelSet;
            obj.nSubCellNodes    = 3;
            obj.nSubCellsByElem  = 3;
            obj.nSubCells        = obj.nElem*obj.nSubCellsByElem;
            obj.nCases           = size(obj.elemCases,2);
        end
        
        function compute(obj)
            obj.computeAllNodesInElem();
            obj.computeNodesInAllSubCells();
            obj.computeIsSubCellInterior();
            obj.computeConnec();
        end
        
        function computeAllNodesInElem(obj)
            vertexNodes = obj.vertexNodesInElem;
            cutNodes    = obj.cutNodesInElem;
            allNodes    = [vertexNodes,cutNodes];
            obj.allNodesInElem = allNodes;
        end
        
        function computeNodesInAllSubCells(obj)
            s.nodesInElem = obj.allNodesInElem;
            s.elemCases = obj.elemCases;              
            s.nSubCellNodes = obj.nSubCellNodes;
            s.nSubCellsByElem = obj.nSubCellsByElem;
            s.nElem = obj.nElem;
            s.nCases = obj.nCases;    
            s.coord = obj.coord;
            s.nSubCells = obj.nSubCells;            
            
            tComputer = TriangleSubCellNodesComputer(s);
            tComputer.compute();
            
            
             obj.allSubCellsConnec = tComputer.allSubCellsConnec;
             obj.isoNode = tComputer.isoNode;
        end
        
   
        function isSubCellInterior = computeIsSubCellInterior(obj)
            
            isoNodeIsFull = obj.levelSet(obj.isoNode) < 0;
            
            isTriangleInterior = isoNodeIsFull;
            
            
            isSubCellInterior = false(obj.nSubCellsByElem,obj.nElem);            
            isSubCellInterior(1,isTriangleInterior)  = true;
            isSubCellInterior(2,~isTriangleInterior) = true;
            isSubCellInterior(3,~isTriangleInterior) = true;
            
            obj.isSubCellInterior = isSubCellInterior;
        end
        
        function computeConnec(obj)
            isInterior = obj.isSubCellInterior(:);
            obj.connec =  obj.allSubCellsConnec(isInterior,:);           
        end
        
        
        
    end
    
end
