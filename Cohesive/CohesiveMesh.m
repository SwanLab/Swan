classdef CohesiveMesh < handle
    
    properties (Access = public)
        fullMesh
        mesh
        centerLineMesh

        listNodeCohesive
        listElemNextCohesive
        listCohesiveElems
        listEdgeCohesive

        pairsMatrix
    end
    
    properties (Access = private)
        baseMesh 
    end
    
    properties (Access = private)
        isNodeCohesive
        isElemCohesive

        separation

        nNodeCohesive
    end
    
    methods (Access = public)
        
        function obj = CohesiveMesh(cParams)
            obj.init(cParams)
            edgesInCohElem = obj.detectFracturedEdges(cParams);
            newCoord    = obj.duplicateNodes();
            centerElemsInEdge         = obj.computeCenterElements();
            normals                           = obj.computeNormals();
            [isLeft, isRight]                 = obj.computeIsLeftIsRight(centerElemsInEdge,normals,edgesInCohElem);
            newCoord    = obj.shiftCoordOfLeftAndRightElements(newCoord,normals);
            newConnec   = obj.updateConnecOfLeftElements(isLeft, newCoord);
            newConnec   = obj.fixFinalElement(cParams,newConnec);

            obj.newMesh(newConnec, newCoord);
            obj.createLineMesh();
        end

        function plot(obj) 
            subplot(1,2,1)
            obj.baseMesh.plot;
            title('BaseMesh')
            subplot(1,2,2)
            obj.fullMesh.plot;
            title('CohesiveMesh')
        end
        
        function createCenterLineMesh(obj)
            coord = obj.mesh.coord;
            nNodes = size(coord,1);
            centerCoord = (coord(1:nNodes-1,:) + coord(2:nNodes,:))/2;
            centerConnec = [(1:nNodes-2)' (2:nNodes-1)'];
            s.coord  = centerCoord;
            s.connec = centerConnec;
            s.kFace  = -1;
            obj.centerLineMesh = Mesh.create(s);
        end

        function A = createMappingMatrix(obj)
            listAllNodesInCoh = [obj.pairsMatrix(:,1);obj.pairsMatrix(:,2)];
            nd = obj.fullMesh.ndim;
            dofsCohGlobal = reshape([nd*listAllNodesInCoh-1; nd*listAllNodesInCoh],1,[]);
            nDofCoh   = length(dofsCohGlobal);
            nDofTotal = obj.fullMesh.nnodes * nd;
            A = sparse(1:nDofCoh, dofsCohGlobal, 1, nDofCoh, nDofTotal);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.separation = 0.0001;
            obj.baseMesh = cParams.baseMesh;
        end

        function edgesInCohElem = detectFracturedEdges(obj, cParams)
            centerEdges      = obj.computeCenterEdge();

            isEdgeCohesive   = cParams.isFracturedLine(centerEdges) & cParams.isFracturedUntil(centerEdges);
            % isEdgeCohesive   = cParams.isFractured(centerEdges);

            obj.listEdgeCohesive = find(isEdgeCohesive);
            nodesInEdgesCohesive = obj.baseMesh.edges.nodesInEdges(obj.listEdgeCohesive,:);
            obj.listNodeCohesive = unique(nodesInEdgesCohesive');

            obj.nNodeCohesive    = length(obj.listNodeCohesive);
            obj.isNodeCohesive   = false(size(obj.baseMesh.coord,1),1);
            obj.isNodeCohesive(obj.listNodeCohesive) = true;

            edgesInElem = obj.baseMesh.edges.edgesInElem;

            obj.isElemCohesive       = any(ismember(edgesInElem,obj.listEdgeCohesive),2);
            obj.listElemNextCohesive = find(obj.isElemCohesive);

            edgesInCohElem = edgesInElem(obj.isElemCohesive,:);
        end

        function centerElemsInCohesiveEdge = computeCenterElements(obj)
            bariCenters = obj.baseMesh.computeBaricenter';
            centerElemsInCohesiveEdge = bariCenters(obj.listElemNextCohesive,:); %Ordenats segons obj.listElemCohesive
        end
        
        function n = computeNormals(obj)
            nodesInEdges = obj.baseMesh.edges.nodesInEdges;
            nodes = nodesInEdges(obj.listEdgeCohesive,:);   
            coords1 = obj.baseMesh.coord(nodes(:,1),:);
            coords2 = obj.baseMesh.coord(nodes(:,2),:);
            swap = (coords2(:,1) < coords1(:,1)) | ...
                    (coords2(:,1) == coords1(:,1) & coords2(:,2) < coords1(:,2));
            tmp = coords1(swap,:);
            coords1(swap,:) = coords2(swap,:);
            coords2(swap,:) = tmp;
            
            t = coords2 - coords1;
            n = [-t(:,2), t(:,1)];
        end

        function [isLeft, isRight] = computeIsLeftIsRight(obj,centerElemsInCohesiveEdge,normals,edgesInCohElem)
            centerEdges  =obj.computeCenterEdge;
            temp             = ismember(edgesInCohElem, obj.listEdgeCohesive);
            cohElemToEdge    = sum(edgesInCohElem .* temp, 2);    % (nElemCoh x 1) 
            centerElem       = centerElemsInCohesiveEdge;         % (nElemCoh x ndim)
            centerEdge       = centerEdges(cohElemToEdge,:);      % (nElemCoh x ndim)
            vectorEdgeToElem = centerElem - centerEdge;

            temp = zeros(size(obj.baseMesh.edges.nodesInEdges,1),2);
            temp(obj.listEdgeCohesive,:) = normals;
            normals = temp(cohElemToEdge,:);

            dotProduct  = sum(vectorEdgeToElem.*normals,2);
            signs = sign(dotProduct);
            signs = signs > 0;
            isLeft = logical(signs);
            isRight = not(isLeft);
        end

        function newCoord = duplicateNodes(obj)
            newCoord    = obj.baseMesh.coord;
            duplicated  = newCoord(obj.isNodeCohesive, :);
            newCoord    = [newCoord; duplicated];
            obj.pairsMatrix = [sort(obj.listNodeCohesive) , linspace(obj.baseMesh.nnodes+1, ...
                obj.baseMesh.nnodes+obj.nNodeCohesive,obj.nNodeCohesive)'];
            obj.isNodeCohesive = [obj.isNodeCohesive; ones(obj.nNodeCohesive,1)];
        end

        function newConnec = updateConnecOfLeftElements(obj,isLeft,newCoord)
            listLeftElems  = obj.listElemNextCohesive(isLeft);
            connec         = obj.baseMesh.connec;
            cohesiveConnec = [obj.pairsMatrix(1:end-1,1), obj.pairsMatrix(2:end,1),   obj.pairsMatrix(2:end,2), obj.pairsMatrix(1:end-1,2)];
            connec         = [connec; cohesiveConnec];

            obj.listCohesiveElems = ((size(connec,1)-size(cohesiveConnec,1)+1):size(connec,1))';
            oldLeftConnec = connec(listLeftElems,:);
            newLeftConnec      = oldLeftConnec;

            idx = ismember(oldLeftConnec,obj.pairsMatrix(:,1));
            newLeftConnec(idx) = arrayfun(@(x) obj.getPair(x), oldLeftConnec(idx));

            connec(listLeftElems,:) = newLeftConnec;
            newConnec = connec;
        end

        function pair = getPair(obj,n)
            pair = obj.pairsMatrix(ismember(obj.pairsMatrix(:,1), n), 2 );
        end
    
        function newMesh(obj,newConnec, newCoord)
            s.connec = newConnec;
            s.coord  = newCoord;
            obj.fullMesh = Mesh.create(s);
        end
    
        function centerEdges = computeCenterEdge(obj)
            obj.baseMesh.computeEdges;
            nodes   = obj.baseMesh.edges.nodesInEdges;    
            coord1  = obj.baseMesh.coord(nodes(:,1),:);
            coord2  = obj.baseMesh.coord(nodes(:,2),:);
            centerEdges = 0.5*(coord1 + coord2);
        end
        
        function newCoord = shiftCoordOfLeftAndRightElements(obj,newCoord,normals)
            shiftVector = 0.5*([normals;zeros(1,size(normals,2))]+ ...
                [zeros(1,size(normals,2));normals]);
            shiftVector(1,:) = normals(1,:); shiftVector(end,:) = normals(end,:);
            shiftVector = shiftVector./vecnorm(shiftVector,2,2); %esta ordenat segons listCohesiveNodes
            
            rightNodes   = obj.listNodeCohesive;
            leftNodes    = getPair(obj,rightNodes);
            newCoord(rightNodes,:) = newCoord(rightNodes,:) - shiftVector*obj.separation/2;
            newCoord(leftNodes,:)  = newCoord(leftNodes,:) + shiftVector*obj.separation/2;   
        end

        function createLineMesh(obj)
            coord     = obj.fullMesh.coord;
            nMidNodes = size(obj.pairsMatrix,1);
            midCoord= (coord(obj.listNodeCohesive,:)+ ...
                coord(obj.pairsMatrix(:,2),:))/2;      
            midConnec =  [(1:nMidNodes-1)' (2:nMidNodes)'];
            s.connec = midConnec;
            s.coord  = midCoord;
            s.kFace  = -1;
            obj.mesh = Mesh.create(s);
        end


        function  newConnec = fixFinalElement(obj,cParams,newConnec)
            centerEdges      = obj.computeCenterEdge();
            listEdgeCohesiveLine = find(cParams.isFracturedLine(centerEdges));
            difference = listEdgeCohesiveLine(not(ismember(listEdgeCohesiveLine,obj.listEdgeCohesive)));
            uniqueEdge = difference(1);
            nodes = obj.baseMesh.edges.nodesInEdges(uniqueEdge,:);
            changedNode = obj.listNodeCohesive(ismember(obj.listNodeCohesive,nodes));
            coord1 = obj.baseMesh.coord(changedNode,:);
            coord2 = obj.baseMesh.coord(nodes(not(ismember(nodes,changedNode))),:);
            t = abs(coord2 - coord1); 
            t = t/norm(t);
            n = [-t(2),t(1)];
            uniqueElems               = find(any(ismember(obj.baseMesh.edges.edgesInElem,uniqueEdge),2));
            bariCenters               = obj.baseMesh.computeBaricenter';
            centerElemsInCohesiveEdge = bariCenters(uniqueElems,:);
            vectorEdgeToElem = centerElemsInCohesiveEdge - centerEdges(uniqueEdge,:);
            signs  = sign(sum(vectorEdgeToElem.*n,2))> 0;
            signs  = signs > 0;
            isLeft = logical(signs);
            uniqueElem = uniqueElems(isLeft);
            changedConnectivity = newConnec(uniqueElem,:);
            changedPosition = ismember(changedConnectivity,changedNode);
            newNode = obj.getPair(changedNode);
            newConnec(uniqueElem,changedPosition) = newNode;
        end
    end
end
