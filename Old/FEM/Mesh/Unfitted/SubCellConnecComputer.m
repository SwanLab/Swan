classdef SubCellConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        connec
    end
    
    properties (Access = private)
        coord
        nElem
        isEdgeCutInElem
        levelSet
        vertexNodesInElem
        
        allNodesInElem
        
        firstCutEdge
        finalNode
        cutNodesInElem
                
        isSubCellInterior
        triangleSubMesher
    end
    
    methods (Access = public)
        
        function obj = SubCellConnecComputer(cParams)
            obj.init(cParams);
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.firstCutEdge      = cParams.firstCutEdge;
            obj.vertexNodesInElem = cParams.backgroundConnec;
            obj.nElem             = cParams.nElem;
            obj.isEdgeCutInElem   = cParams.isEdgeCutInElem;
            obj.coord             = cParams.coord;
            obj.levelSet          = cParams.levelSet;
            obj.finalNode         = cParams.finalNode;
            
        end
        
        function compute(obj)
            obj.computeCutNodePerElem();
            obj.computeAllNodesInElem();
            obj.computeNodesInAllSubCells();
            obj.computeIsSubCellsInterior();
            obj.computeConnec();
        end
        
        function computeAllNodesInElem(obj)
            vertexNodes = obj.vertexNodesInElem;
            cutNodes    = obj.cutNodesInElem;
            allNodes    = [vertexNodes,cutNodes];
            obj.allNodesInElem = allNodes;
        end
        
        function computeNodesInAllSubCells(obj)
            s.nodesInElem     = obj.allNodesInElem;
            s.isEdgeCutInElem = obj.isEdgeCutInElem;
            s.nElem           = obj.nElem;
            s.coord           = obj.coord;
            tComputer = TriangleSubCellNodesComputer(s);
            tComputer.compute();
            obj.triangleSubMesher = tComputer;
        end
           
        function computeIsSubCellsInterior(obj)
            isTriInt  = obj.computeIsSubCellTriangleInterior();
            nSubCells = obj.triangleSubMesher.nSubCellsByElem;
            itIs = false(nSubCells,obj.nElem);
            itIs(1,isTriInt)  = true;
            itIs(2,~isTriInt) = true;
            itIs(3,~isTriInt) = true;
            obj.isSubCellInterior = itIs;
        end
        
        function itIs = computeIsSubCellTriangleInterior(obj)
            isoNode = obj.triangleSubMesher.isoNode;
            isoNodeIsFull = obj.levelSet(isoNode) < 0;
            itIs = isoNodeIsFull;
        end
        
        function computeConnec(obj)
            allConnec  = obj.triangleSubMesher.allSubCellsConnec;
            isInterior = obj.isSubCellInterior(:);
            obj.connec = allConnec(isInterior,:);
        end
        
        function cutNodePerElemen = computeCutNodePerElem(obj)
            firstCutEdgePerElem = obj.firstCutEdge;
            cutNodePerElemen = firstCutEdgePerElem + obj.finalNode;
            obj.cutNodesInElem = cutNodePerElemen;
        end
        
        
    end
    
end
