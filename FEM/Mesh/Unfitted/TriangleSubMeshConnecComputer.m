classdef TriangleSubMeshConnecComputer < handle
    
    properties (Access = public)
        connecSubTriangle
        connecSubQuad
        nSubCellNodes
        nSubCellsByElem
    end
    
    properties (Access = private)
        coord
        
        localTriangleNodeCases
        localQuadNodeCases
        
        cellNodes
        nSubCellsByQuad
        nSubCases
        nElemInCase       
    end
    
    methods (Access = public)
        
        function obj = TriangleSubMeshConnecComputer(cParams)
            obj.init(cParams)
        end
        
        function nodesInSubCells = compute(obj,nodes,icase)
            obj.nElemInCase = size(nodes,1);
            obj.cellNodes   = nodes;
            nodesInSubCells = zeros(obj.nSubCellsByElem,obj.nSubCellNodes,obj.nElemInCase); 
            obj.connecSubTriangle = obj.computeSubTriangleConnec(icase);
            obj.connecSubQuad     = obj.computeSubQuadConnec(icase);
            
            nodesInSubCells(1,:,:)   = obj.connecSubTriangle;
            nodesInSubCells(2:3,:,:) = obj.connecSubQuad;
        end
        
       
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord = cParams.coord;            
            obj.nSubCellNodes   = 3;
            obj.nSubCellsByQuad = 2;
            obj.nSubCases       = 2;
            obj.nSubCellsByElem = 3;
            obj.computeTriangleLocalNodeCases();
            obj.computeLocalNodesCases();
        end
        
        function computeTriangleLocalNodeCases(obj)
            nodes = [1 4 5;4 2 5; 4 3 5];
            obj.localTriangleNodeCases = nodes;
        end
        
        function computeLocalNodesCases(obj)
            nodes(:,:,1,1) = [4 2 3;5 4 3];
            nodes(:,:,2,1) = [2 3 5;4 2 5];
            
            nodes(:,:,1,2) = [4 5 3;1 4 3];
            nodes(:,:,2,2) = [1 4 5;1 5 3];
            
            nodes(:,:,1,3) = [1 2 5;2 4 5];
            nodes(:,:,2,3) = [1 4 5;1 2 4];            
            obj.localQuadNodeCases = nodes;
        end
        
        function nodesT = computeSubTriangleConnec(obj,icase)
            localNodes   = obj.localTriangleNodeCases(icase,:);
            nodesT = zeros(obj.nSubCellNodes,obj.nElemInCase);
            for inode = 1:obj.nSubCellNodes
                localNode = localNodes(inode);
                node = obj.cellNodes(:,localNode);
                nodesT(inode,:) = node;
            end
        end
        
        function nodeQ = computeSubQuadConnec(obj,icase)
            localNodes = obj.localQuadNodeCases(:,:,:,icase);
            nodesSubCases = obj.computeNodesSubCases(localNodes);
            nodeQ = obj.computeBestCase(nodesSubCases);
        end        
        
        function nodeQ = computeBestCase(obj,nodesSubCases)
           s.coord         = obj.coord;
           s.nodesSubCases = nodesSubCases;
           bestCase = BestSubCellCaseSelector(s);
           nodeQ = bestCase.compute();
        end
        
        function nodesSubCases = computeNodesSubCases(obj,localNodesCase)
            nodesSubCases = zeros(obj.nSubCellsByQuad,obj.nSubCellNodes,obj.nSubCases,obj.nElemInCase);
            for isubCase = 1:obj.nSubCases
                localNodes = localNodesCase(:,:,isubCase);
                nodesT = zeros(obj.nSubCellsByQuad,obj.nSubCellNodes,obj.nElemInCase);
                for isubCell = 1:obj.nSubCellsByQuad
                    for inode = 1:obj.nSubCellNodes
                        localNode = localNodes(isubCell,inode);
                        node = obj.cellNodes(:,localNode);
                        nodesT(isubCell,inode,:) = node;
                    end
                end
                nodesSubCases(:,:,isubCase,:) = nodesT;
            end
        end
        
    
        
    end
    
    
    
end