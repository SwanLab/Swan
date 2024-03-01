 classdef AllNodesInElemCoordinatesComputer < handle
    
     properties (Access = public)
        xAllNodesInElem
        xCutInElem
        nodesInCutEdges
     end
     
    properties (Access = private)
        
    end
     
    properties (Access = private)
        interpolation
        xIsoNodes
        nElem
        nEdgeByElem
        localNodeByEdgeByElem
        xCutEdgePoint
        edgeCutPointInElem
        all2Cut
    end
    
    methods (Access = public)
       
        function obj = AllNodesInElemCoordinatesComputer(cParams)
            obj.init(cParams);
            obj.computeXisoNodes();
        end
        
        function compute(obj)
            obj.initXnodesElem();
            obj.addXisoInXall();
            obj.computeNodesInCutEdges()
            obj.computeXcut(); 
            obj.addXcutInXall();
        end
        
    end
     
    methods (Access = private)
        
        function init(obj,cParams)
       %     obj.nCutEdgeByElem        = cParams.nCutEdgeByElem;
            obj.localNodeByEdgeByElem = cParams.localNodeByEdgeByElem;
            obj.xCutEdgePoint         = cParams.xCutEdgePoint;
            obj.edgeCutPointInElem    = cParams.edgeCutPointInElem;
            obj.all2Cut               = cParams.all2Cut;
            obj.nElem                 = size(obj.localNodeByEdgeByElem,1);
            obj.nEdgeByElem           = size(obj.localNodeByEdgeByElem,2);
            obj.createInterpolation();
        end
        
        function createInterpolation(obj)
            type = 'LINE';
            int = Interpolation.create(type,'LINEAR');
            obj.interpolation = int;
        end
        
        function computeXisoNodes(obj)
            switch obj.nEdgeByElem
                case 3
                    type = 'TRIANGLE';
                case 6
                    type = 'TETRAHEDRA';
            end
            int = Interpolation.create(type,'LINEAR');
            obj.xIsoNodes = int.pos_nodes;
        end
        
        function addXisoInXall(obj)
            nNode = size(obj.xIsoNodes,1);
            xIso = repmat(obj.xIsoNodes',[1 1 obj.nElem]);
            for inode = 1:nNode
                obj.xAllNodesInElem(:,inode,:) = xIso(:,inode,:);
            end
        end
        
        function addXcutInXall(obj)
            nNode = size(obj.xIsoNodes,1);  
            xC = obj.xCutInElem;
            for inode = 1:obj.all2Cut.nCutEdgeByElem
                obj.xAllNodesInElem(:,inode+nNode,:) = xC(:,inode,:);
            end
        end
        
        function initXnodesElem(obj)
            nNode = size(obj.xIsoNodes,1);
            nDim  = size(obj.xIsoNodes,2);
            nAllNodes = nNode;
            nodes = zeros(nDim,nAllNodes,obj.nElem);
            obj.xAllNodesInElem = nodes;
        end
        
        function computeNodesInCutEdges(obj)
            allNodesByEdge = obj.localNodeByEdgeByElem;
            nNodeByEdge = size(allNodesByEdge,3);
            nodesCutEdges = zeros(obj.nElem,obj.all2Cut.nCutEdgeByElem,nNodeByEdge);
            for inode = 1:nNodeByEdge
                nodes = allNodesByEdge(:,:,inode);
                nodesCutEdges(:,:,inode) = obj.all2Cut.compute(nodes);
            end
            obj.nodesInCutEdges = nodesCutEdges;
        end
        
        function computeXcut(obj)
            xC = obj.initXcut();
            for iedge = 1:obj.all2Cut.nCutEdgeByElem
                [xA,xB] = obj.computeXnodes(iedge);
                [shapeA,shapeB] = obj.computeShapes(iedge);
                x = xA.*shapeA + xB.*shapeB; 
                xC(:,iedge,:) = x;
            end
            obj.xCutInElem = xC;
        end
        
        function xCut = initXcut(obj)
            nDim = size(obj.xIsoNodes,2);
            xCut = zeros(nDim,obj.all2Cut.nCutEdgeByElem,obj.nElem);
        end
        
        function [xA,xB] = computeXnodes(obj,iedge)
            nodes = obj.nodesInCutEdges;
            nodeA = nodes(:,iedge,1);
            nodeB = nodes(:,iedge,2);
            xA = obj.xIsoNodes(nodeA,:)';
            xB = obj.xIsoNodes(nodeB,:)';
        end
        
        function [shapeA,shapeB] = computeShapes(obj,iedge)
            edge = obj.edgeCutPointInElem(:,iedge);
            xCutPoint(1,:)  = obj.xCutEdgePoint(edge,1);
            shapes = obj.interpolation.computeShapeFunctions(xCutPoint);
            shapeA = shapes(1,:);
            shapeB = shapes(2,:);
        end

    end
     
 end
