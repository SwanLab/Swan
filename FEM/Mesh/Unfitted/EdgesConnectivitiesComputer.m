classdef EdgesConnectivitiesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        nodesInEdges
        edgesInElem
        nEdgeByElem
        nNodeByEdge   
        localNodeByEdgeByElem
    end
    
    properties (Access = private)
        nodesByElem
        nElem
        nAllEdges
        
        localEdgesInElem
        allToUnique
        uniqueToAll
        nodesInAllEdges
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
            obj.computeLocalNodeByEdgeByElem();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodesByElem = cParams.nodesByElem;
            obj.localEdgesInElem =  [1 2; 2 3; 3 1];
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
        
        function computeLocalNodeByEdgeByElem(obj)
            for iedge = 1:3
            edges = obj.edgesInElem(:,iedge);
            node(:,1) = obj.nodesInEdges(edges,1);
            node(:,2) = obj.nodesInEdges(edges,2);


            localNodeA = obj.localEdgesInElem(iedge,1);
            localNodeB = obj.localEdgesInElem(iedge,2);

            node2(:,1) = obj.nodesByElem(:,localNodeA);
            node2(:,2) = obj.nodesByElem(:,localNodeB);
            
            isLocalAndGlobalEquallyOriented = sum(abs(node - node2),2) < 1;
            itIs = isLocalAndGlobalEquallyOriented;
            localNodeByEdgeByElem(itIs,1,iedge)  = localNodeA;
            localNodeByEdgeByElem(itIs,2,iedge)  = localNodeB;
            localNodeByEdgeByElem(~itIs,1,iedge) = localNodeB;
            localNodeByEdgeByElem(~itIs,2,iedge) = localNodeA;
            end
            
            obj.localNodeByEdgeByElem = localNodeByEdgeByElem;
            
        end
        
    end
end