classdef FacesConnectivitiesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        nodesInFaces
        facesInElem
        nFaceByElem
        nNodeByFace
        localNodeByFaceByElem
    end
    
    properties (Access = private)
        allToUnique
        uniqueToAll
        nodesInAllFaces
    end
    
    properties (Access = private)
        nodesByElem
        nElem
        nAllFaces
        localFacesInElem
        type
    end
    
    methods (Access = public)
        
        function obj = FacesConnectivitiesComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeNodesInAllFaces();
            obj.computeUniqueFacesIndex();
            obj.computeNodesInUniqueFaces();
            obj.computeFacesInElem();
            obj.computeLocalOrientedFaceConnec();
        end
        
        function cV = computeConnectedVertex(obj,vertex)
            e  = obj.computeAllEdgesOfVertex(vertex);
            cV = obj.computeOtherVertexOfEdge(e,vertex);
        end         
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodesByElem = cParams.nodesByElem;
            obj.type = cParams.type;
            switch obj.type
                case 'TRIANGLE'
                    obj.localFacesInElem = [];
                    
                case 'QUAD'
                    obj.localFacesInElem = [];

                case 'TETRAHEDRA'
                    obj.localFacesInElem = [1 2 4; 1 3 4; 1 2 3; 2 3 4];
                    
                case 'LINE'
                    obj.localFacesInElem =  [];
                    
                case 'HEXAHEDRA'
                    obj.localFaceInElem  = [1 2; 4 1; 1 5; 2 3; 2 6; 3 4; ...
                        3 7; 4 8; 5 6; 8 5; 6 7; 7 8];
            end
            obj.nElem = size(obj.nodesByElem,1);
            obj.nFaceByElem = size(obj.localFacesInElem,1);
            obj.nNodeByFace = size(obj.localFacesInElem,2);
            obj.nAllFaces = obj.nElem*obj.nFaceByElem;
        end
        
        function computeNodesInAllFaces(obj)
            nodesInF = zeros(obj.nAllFaces,obj.nNodeByFace);
            allElem = 1:obj.nElem;
            for inode = 1:obj.nNodeByFace
                for iedge = 1:obj.nFaceByElem
                    lNode = obj.localFacesInElem(iedge,inode);
                    index = (allElem - 1)*obj.nFaceByElem + iedge;
                    nodes = obj.nodesByElem(:,lNode);
                    nodesInF(index,inode) = nodes;
                end
            end
            obj.nodesInAllFaces = nodesInF;
        end
        
        function computeUniqueFacesIndex(obj)
            nE = obj.nodesInAllFaces;
            [nE] = sort(nE,2);
            [~,a2u,u2a] = unique(nE,'rows');
            obj.allToUnique = a2u;
            obj.uniqueToAll = u2a;
        end
        
        function computeNodesInUniqueFaces(obj)
            a2u = obj.allToUnique;
            nodesInUniqueFaces = obj.nodesInAllFaces(a2u,:);
            obj.nodesInFaces = nodesInUniqueFaces;
        end
        
        function computeFacesInElem(obj)
            uniqueFaceByElem = obj.computeUniqueFacesByElem();
            fInElem = zeros(obj.nElem,obj.nFaceByElem);
            for iFace = 1:obj.nFaceByElem
                face = uniqueFaceByElem(:,iFace);
                fInElem(:,iFace) = obj.uniqueToAll(face);
            end
            obj.facesInElem = fInElem;
        end
        
        function edgeByElem = computeUniqueFacesByElem(obj)
            edges(:,1) = 1:obj.nAllFaces;
            edgeByElem = reshape(edges,obj.nFaceByElem,obj.nElem);
            edgeByElem = edgeByElem';
        end
        
        function computeLocalOrientedFaceConnec(obj)
            faceConnec = zeros(obj.nElem,obj.nFaceByElem,obj.nNodeByFace);
            for iFace = 1:obj.nFaceByElem
                isOriented = obj.isFaceOriented(iFace);
                nodeA = obj.localFacesInElem(iFace,1); % CHANGE THIS
                nodeB = obj.localFacesInElem(iFace,2);
                faceConnec(isOriented,iFace,1)  = nodeA;
                faceConnec(isOriented,iFace,2)  = nodeB;
                faceConnec(~isOriented,iFace,1) = nodeB;
                faceConnec(~isOriented,iFace,2) = nodeA;
            end
            obj.localNodeByFaceByElem = faceConnec;
        end
        
        function itIs = isFaceOriented(obj,iFace)
            globalNodes = obj.computeGlobalNode(iFace);
            localNodes  = obj.computeLocalNodes(iFace);
            itIs = obj.isEquallyOriented(globalNodes,localNodes);
        end
        
        function globalNodes = computeGlobalNode(obj,iFace)
            faces       = obj.facesInElem(:,iFace);
            globalNodes = obj.nodesInFaces(faces,:);
        end
        
        function localNodes = computeLocalNodes(obj,iFace)
            localNode   = obj.localFacesInElem(iFace,:);
            localNodes  = obj.nodesByElem(:,localNode);
        end
        
        function edges = computeAllEdgesOfVertex(obj,vertex)
            vertexInEdges = obj.nodesInEdges;            
            isInEdge = vertexInEdges == vertex;            
            allEdges  = 1:size(isInEdge,1);
            allEdgesA(:,1) = allEdges(isInEdge(:,1));
            allEdgesB(:,1) = allEdges(isInEdge(:,2));
            edges = [allEdgesA;allEdgesB];            
        end        
        
        function oV = computeOtherVertexOfEdge(obj,edge,vertex)
            vertexesInEdges = obj.nodesInEdges;
            vertexesOfEdge  = vertexesInEdges(edge,:);
            oVB = setdiff(vertexesOfEdge(:,2),vertex);
            oVA = setdiff(vertexesOfEdge(:,1),vertex);
            oV = [oVB;oVA];
        end            
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isEquallyOriented(gNode,lNode)
            itIs = sum(abs(gNode - lNode),2) < 1;
        end
        
    end
end