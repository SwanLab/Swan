classdef EdgesConnectivitiesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        nodesInEdges
        edgesInElem
        nEdgeByElem
        nNodeByEdge
        localNodeByEdgeByElem
    end
    
    properties (Access = private)
        allToUnique
        uniqueToAll
        nodesInAllEdges                
    end
    
    properties (Access = private)
        nodesByElem
        nElem
        nAllEdges       
        localEdgesInElem
    end
    
    methods (Access = public)
        
        function obj = EdgesConnectivitiesComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeNodesInAllEdges();
            obj.computeUniqueEdgesIndex();
            obj.computeNodesInUniqueEdges();
            obj.computeEdgesInElem();
            obj.computeLocalOrientedEdgeConnec();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodesByElem = cParams.nodesByElem;
            nNodes = size(cParams.nodesByElem,2);
            switch nNodes
                case 3            
                    obj.localEdgesInElem = nchoosek(1:nNodes,2);%[1 2; 2 3; 3 1];
                    obj.localEdgesInElem = [1 2; 2 3; 3 1];                    
                case 4
                    obj.localEdgesInElem =  nchoosek(1:nNodes,2);%[1 2; 2 3; 3 1];
            end            
            obj.nElem = size(obj.nodesByElem,1);
            obj.nEdgeByElem = size(obj.localEdgesInElem,1);
            obj.nNodeByEdge = size(obj.localEdgesInElem,2);
            obj.nAllEdges = obj.nElem*obj.nEdgeByElem;
        end
        
        function computeNodesInAllEdges(obj)
            nodesInE = zeros(obj.nAllEdges,obj.nNodeByEdge);
            allElem = 1:obj.nElem;
            for inode = 1:obj.nNodeByEdge
                for iedge = 1:obj.nEdgeByElem
                    lNode = obj.localEdgesInElem(iedge,inode);
                    index = (allElem - 1)*obj.nEdgeByElem + iedge;
                    nodes = obj.nodesByElem(:,lNode);
                    nodesInE(index,inode) = nodes;
                end
            end
            obj.nodesInAllEdges = nodesInE;
        end
        
        function computeUniqueEdgesIndex(obj)
            nE = obj.nodesInAllEdges;
            [nE] = sort(nE,2);
            [~,a2u,u2a] = unique(nE,'rows');
            obj.allToUnique = a2u;
            obj.uniqueToAll = u2a;
        end
        
        function computeNodesInUniqueEdges(obj)
            a2u = obj.allToUnique;
            nodesInUniqueEdges = obj.nodesInAllEdges(a2u,:);
            obj.nodesInEdges = nodesInUniqueEdges;
        end
        
        function computeEdgesInElem(obj)
            uniqueEdgeByElem = obj.computeUniqueEdgesByElem();
            eInElem = zeros(obj.nElem,obj.nEdgeByElem);
            for iEdge = 1:obj.nEdgeByElem
                edge = uniqueEdgeByElem(:,iEdge);
                eInElem(:,iEdge) = obj.uniqueToAll(edge);
            end
            obj.edgesInElem = eInElem;
        end
        
        function edgeByElem = computeUniqueEdgesByElem(obj)
            edges(:,1) = 1:obj.nAllEdges;
            edgeByElem = reshape(edges,obj.nEdgeByElem,obj.nElem);
            edgeByElem = edgeByElem';
        end
        
        function computeLocalOrientedEdgeConnec(obj)
            edgeConnec = zeros(obj.nElem,obj.nEdgeByElem,obj.nNodeByEdge);
            for iedge = 1:obj.nEdgeByElem
                isOriented = obj.isEdgeOriented(iedge);                
                nodeA = obj.localEdgesInElem(iedge,1);
                nodeB = obj.localEdgesInElem(iedge,2);
                edgeConnec(isOriented,iedge,1)  = nodeA;
                edgeConnec(isOriented,iedge,2)  = nodeB;
                edgeConnec(~isOriented,iedge,1) = nodeB;
                edgeConnec(~isOriented,iedge,2) = nodeA;
            end            
            obj.localNodeByEdgeByElem = edgeConnec;            
        end
        
        function itIs = isEdgeOriented(obj,iedge)
            globalNodes = obj.computeGlobalNode(iedge);
            localNodes  = obj.computeLocalNodes(iedge);
            itIs = obj.isEquallyOriented(globalNodes,localNodes);
        end
        
        function globalNodes = computeGlobalNode(obj,iedge)
            edges       = obj.edgesInElem(:,iedge);
            globalNodes = obj.nodesInEdges(edges,:);
        end
        
        function localNodes = computeLocalNodes(obj,iedge)
            localNode   = obj.localEdgesInElem(iedge,:);
            localNodes  = obj.nodesByElem(:,localNode);
        end
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isEquallyOriented(gNode,lNode)
            itIs = sum(abs(gNode - lNode),2) < 1;
        end
        
    end
end