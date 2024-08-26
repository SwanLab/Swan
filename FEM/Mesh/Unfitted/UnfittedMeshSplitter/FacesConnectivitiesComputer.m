classdef FacesConnectivitiesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        nodesInFaces
        facesInElem
        nFaceByElem
        nNodeByFace
        nEdgeByFace
%         localNodeByFaceByElem
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
%             obj.computeLocalOrientedFaceConnec();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodesByElem = cParams.nodesByElem;
            obj.type = cParams.type;
            switch obj.type
                case 'TRIANGLE'
                    obj.localFacesInElem = [1 2 3];
                    obj.nEdgeByFace = 3;
                    
                case 'QUAD'
                    obj.localFacesInElem = [1 2 3 4];
                    obj.nEdgeByFace = 4;

                case 'TETRAHEDRA'
                    obj.localFacesInElem = [1 2 4; 1 3 4; 1 2 3; 2 3 4];
                    obj.nEdgeByFace = 3;
                    
                case 'LINE'
                    obj.localFacesInElem =  [1 2];
                    obj.nEdgeByFace = 1;
                    
                case 'HEXAHEDRA'
                    obj.localFacesInElem  = [1 2 3 4; 1 2 6 5; 1 4 8 5; ...
                        2 3 7 6; 3 4 8 7; 5 6 7 8];
                    obj.nEdgeByFace = 4;
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
                nodesRotated = obj.isFaceRotated(iFace);
                
                nodeA = obj.localFacesInElem(iFace,1); % CHANGE THIS
                nodeB = obj.localFacesInElem(iFace,2);
                faceConnec(nodesRotated,iFace,1)  = nodeA;
                faceConnec(nodesRotated,iFace,2)  = nodeB;
                faceConnec(~nodesRotated,iFace,1) = nodeB;
                faceConnec(~nodesRotated,iFace,2) = nodeA;
            end
            obj.localNodeByFaceByElem = faceConnec;
        end
        
        function howMuch = isFaceRotated(obj,iFace)
            globalNodes = obj.computeGlobalNode(iFace);
            localNodes  = obj.computeLocalNodes(iFace);
            howMuch = obj.isEquallyOriented(globalNodes,localNodes);
            if ~all(howMuch==1)
                disp('hola')
            end
        end
        
        function globalNodes = computeGlobalNode(obj,iFace)
            faces       = obj.facesInElem(:,iFace);
            globalNodes = obj.nodesInFaces(faces,:);
        end
        
        function localNodes = computeLocalNodes(obj,iFace)
            localNode   = obj.localFacesInElem(iFace,:);
            localNodes  = obj.nodesByElem(:,localNode);
        end         
        
    end
    
    methods (Access = private, Static)
        
        function howMuch = isEquallyOriented(gNode,lNode)
            itIs = sum(abs(gNode - lNode),2) < 1;
            howMuch = size(itIs);
            for i = 1:length(itIs)
                if itIs(i)==0
                    
                else
                    howMuch(i) = 0;
                end
            end
        end
        
    end
end