classdef FacesConnectivitiesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        nodesInFaces
        facesInElem
        nFaceByElem
        nNodeByFace
        nEdgeByFace
        localNodeByFaceByElem
        localFacesInElem
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

        function sides = getFacesOrientation(obj)
            l_f = obj.localNodeByFaceByElem;
            r_f = obj.localFacesInElem;

            sides = zeros(size(obj.facesInElem));

            for iFace = 1:obj.nFaceByElem
                for iRot = 1:obj.nNodeByFace
                    glo = squeeze(l_f(:,iFace,:));
                    ref = circshift(r_f(iFace,:),iRot-1);
                    sides(all(glo==ref,2),iFace) = 1;
                end
                sides(sides(:,iFace)==0,iFace) = -1;
            end
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
            cases = perms(1:obj.nNodeByFace);

            for iFace = 1:obj.nFaceByElem
                nodes = obj.localFacesInElem(iFace,:);

                for iCase = 1:factorial(obj.nNodeByFace)
                    nodesRotated = obj.isFaceRotated(iFace, iCase);
                    faceConnec(nodesRotated == 1,iFace,:) = repmat(nodes(cases(iCase,:)),[sum(nodesRotated) 1 1]);
                end
            end
            obj.localNodeByFaceByElem = faceConnec;
        end
        
        function howMuch = isFaceRotated(obj,iFace,iCase)
            globalNodes = obj.computeGlobalNode(iFace);
            localNodes  = obj.computeLocalNodes(iFace,iCase);
            howMuch = obj.isEquallyOriented(globalNodes,localNodes);
        end
        
        function globalNodes = computeGlobalNode(obj,iFace)
            faces       = obj.facesInElem(:,iFace);
            globalNodes = obj.nodesInFaces(faces,:);
        end
        
        function localNodes = computeLocalNodes(obj,iFace,iCase)
            cases = perms(1:obj.nNodeByFace);
            localNode   = obj.localFacesInElem(iFace,cases(iCase,:));
            localNodes  = obj.nodesByElem(:,localNode);
        end
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isEquallyOriented(gNode,lNode)
            itIs = sum(abs(gNode - lNode),2) < 1;
        end
        
    end
end