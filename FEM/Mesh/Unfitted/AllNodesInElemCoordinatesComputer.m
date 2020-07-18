 classdef AllNodesInElemCoordinatesComputer < handle
    
     properties (Access = public)
        xAllNodesInElem   
        xCut
     end
     
    properties (Access = private)
        
    end
     
    properties (Access = private)
        interpolation
        xIsoNodes   
        nCutEdgeByElem
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
            obj.computeXcut(); 
            obj.addXcutInXall();
        end
        
    end
     
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nCutEdgeByElem        = cParams.nCutEdgeByElem;
            obj.localNodeByEdgeByElem = cParams.localNodeByEdgeByElem;
            obj.xCutEdgePoint         = cParams.xCutEdgePoint;
            obj.edgeCutPointInElem    = cParams.edgeCutPointInElem;
            obj.all2Cut               = cParams.all2Cut;
            obj.nElem                 = size(obj.localNodeByEdgeByElem,1);
            obj.nEdgeByElem           = size(obj.localNodeByEdgeByElem,2);            
            obj.createInterpolation();
        end
        
        function createInterpolation(obj)
            m.type = 'LINE';
            m.coord  = [];
            m.connec = [];
            int = Interpolation.create(m,'LINEAR');
            obj.interpolation = int;
        end
        
        function computeXisoNodes(obj)
            m.type = 'TRIANGLE';
            m.coord  = [];
            m.connec = [];
            int = Interpolation.create(m,'LINEAR');
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
            xC = obj.xCut;
            for inode = 1:obj.nCutEdgeByElem
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
        
        function computeXcut(obj)
            nodes = obj.computeIsoCutNodesInElem();
            xC = obj.initXcut();
            for iedge = 1:2                
                [xA,xB] = obj.computeXnodes(nodes,iedge);
                [shapeA,shapeB] = obj.computeShapes(iedge);
                x = xA.*shapeA + xB.*shapeB; 
                xC(:,iedge,:) = x;
            end
            obj.xCut = xC;
        end
        
        function nodesCutEdges = computeIsoCutNodesInElem(obj)
            allNodesByEdge = obj.localNodeByEdgeByElem;
            nNodeByEdge = size(allNodesByEdge,3);
            nodesCutEdges = zeros(obj.nElem,obj.nCutEdgeByElem,nNodeByEdge);
            for inode = 1:nNodeByEdge
                nodes = allNodesByEdge(:,:,inode);
                nodesCutEdges(:,:,inode) = obj.all2Cut.compute(nodes);
            end
        end        
        
        function xCut = initXcut(obj)
            nDim = size(obj.xIsoNodes,2);            
            xCut = zeros(nDim,obj.nCutEdgeByElem,obj.nElem);            
        end
        
        function [xA,xB] = computeXnodes(obj,nodes,iedge)
            nodeA = nodes(:,iedge,1);
            nodeB = nodes(:,iedge,2);
            xA = obj.xIsoNodes(nodeA,:)';
            xB = obj.xIsoNodes(nodeB,:)';
        end
        
        function [shapeA,shapeB] = computeShapes(obj,iedge)
            edge = obj.edgeCutPointInElem(:,iedge);
            xCutPoint(1,:)  = obj.xCutEdgePoint(edge,1);          
            obj.interpolation.computeShapeDeriv(xCutPoint);
            shapes = obj.interpolation.shape;
            shapeA = shapes(1,:);
            shapeB = shapes(2,:);
        end       

    end    
     
 end




