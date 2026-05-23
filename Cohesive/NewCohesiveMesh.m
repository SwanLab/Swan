classdef NewCohesiveMesh < handle
    
    properties (Access = public)
        fullMesh
        mesh
        centerLineMesh

        listNodeCohesive

        pairsMatrixReal
        pairsMatrixFull
    end
    
    properties (Access = private)
        baseMesh 

        isEdgeCohesiveReal
        isEdgeCohesiveAuxiliar

        separation
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = NewCohesiveMesh(cParams)
            obj.init(cParams)
            obj.detectFracturedEdges(cParams);
            newCoord    = obj.duplicateNodes();
            centerElemsInEdgeReal         = obj.computeCenterElements();
            [normalsOfEdgesFull, normalsOfEdgeReal]  = obj.computeNormals();
            [isLeft, isRight]                 = obj.computeIsLeftIsRight(centerElemsInEdgeReal,normalsOfEdgesFull);
            newCoord    = obj.shiftCoordOfLeftAndRightElements(newCoord,normalsOfEdgesFull);
            newConnec   = obj.updateConnecOfLeftElements(isLeft, newCoord);
            % newConnec   = obj.fixFinalElement(cParams,newConnec);

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
            listAllNodesInCoh = [obj.pairsMatrixFull(:,1);obj.pairsMatrixFull(:,2)];
            ndim = obj.fullMesh.ndim;
            dofsCohGlobal = reshape([ndim*listAllNodesInCoh-1; ndim*listAllNodesInCoh],1,[]);
            nDofCoh   = length(dofsCohGlobal);
            nDofTotal = obj.fullMesh.nnodes * ndim;
            A = sparse(1:nDofCoh, dofsCohGlobal, 1, nDofCoh, nDofTotal);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.separation = 0.0001;
            obj.baseMesh = cParams.baseMesh;
        end

        function edgesInCohElemFull = detectFracturedEdges(obj, cParams)
            centerEdges      = obj.computeCenterEdge();

            x = cParams.xCohLineMax;
            y = cParams.yCohLine;
            isEdgeCohesiveFull  = abs(centerEdges(:,2) - y) <= 1e-10;
            isFracturedUntil = centerEdges(:,1) - x  <= 1e-10;

            obj.isEdgeCohesiveReal        = isEdgeCohesiveFull & isFracturedUntil;
            obj.isEdgeCohesiveAuxiliar    = isEdgeCohesiveFull & not(isFracturedUntil);

            % edgesInElem = obj.baseMesh.edges.edgesInElem;
            % isElemCohesiveReal  = any(ismember(edgesInElem,obj.listEdgeCohesiveReal),2);
            % edgesInCohElemFull = edgesInElem(obj.isElemCohesiveFull,:);

            
            % obj.listEdgeCohesiveReal = find(isEdgeCohesiveFull);
            % nodesInEdgesCohesive = obj.baseMesh.edges.nodesInEdges(obj.listEdgeCohesiveReal,:);
            % obj.listNodeCohesive = unique(nodesInEdgesCohesive');
            % 
            % obj.nNodeCohesive    = length(obj.listNodeCohesive);
            % obj.isNodeCohesive   = false(size(obj.baseMesh.coord,1),1);
            % obj.isNodeCohesive(obj.listNodeCohesive) = true;
            % 
            % edgesInElem = obj.baseMesh.edges.edgesInElem;
            % 
            % obj.isElemCohesive       = any(ismember(edgesInElem,obj.listEdgeCohesiveReal),2);
            % obj.listElemNextCohesive = find(obj.isElemCohesive);
            % 
            % edgesInCohElem = edgesInElem(obj.isElemCohesive,:);
        end

        function centerElemsInCohesiveEdgeFull = computeCenterElements(obj)
            edgesInElem = obj.baseMesh.edges.edgesInElem;

            listEdgeCohesiveReal = find(obj.isEdgeCohesiveReal);
            isElemCohesiveReal       = any(ismember(edgesInElem,listEdgeCohesiveReal),2);
            listElemNextCohesiveReal = find(isElemCohesiveReal);

            listEdgeCohesiveFull = find(obj.isEdgeCohesiveReal | obj.isEdgeCohesiveAuxiliar);
            isElemCohesiveFull       = any(ismember(edgesInElem,listEdgeCohesiveFull),2);
            listElemNextCohesiveFull = find(isElemCohesiveFull);


            bariCenters = obj.baseMesh.computeBaricenter';
            centerElemsInCohesiveEdgeReal = bariCenters(listElemNextCohesiveReal,:); %Ordenats segons obj.listElemCohesive
            centerElemsInCohesiveEdgeFull = bariCenters(listElemNextCohesiveFull,:);
        end
        
        function [nf, nr] = computeNormals(obj)
            nodesInEdges = obj.baseMesh.edges.nodesInEdges;
            listEdgeCohesiveFull = find(obj.isEdgeCohesiveReal | obj.isEdgeCohesiveAuxiliar);
            listEdgeCohesiveReal = find(obj.isEdgeCohesiveReal);

            nodes = nodesInEdges(listEdgeCohesiveFull,:);   
            coords1 = obj.baseMesh.coord(nodes(:,1),:);
            coords2 = obj.baseMesh.coord(nodes(:,2),:);
            swap = (coords2(:,1) < coords1(:,1)) | ...
                    (coords2(:,1) == coords1(:,1) & coords2(:,2) < coords1(:,2));
            tmp = coords1(swap,:);
            coords1(swap,:) = coords2(swap,:);
            coords2(swap,:) = tmp;
            
            t = coords2 - coords1;
            nf = [-t(:,2), t(:,1)];
            nr = nf(ismember(listEdgeCohesiveFull,listEdgeCohesiveReal),:);
        end

        function [isLeft, isRight] = computeIsLeftIsRight(obj,centerElemsInCohesiveEdge,normalsFull)
            centerEdges      = obj.computeCenterEdge;

            edgesInElem = obj.baseMesh.edges.edgesInElem;
            listEdgeCohesiveFull = find(obj.isEdgeCohesiveReal | obj.isEdgeCohesiveAuxiliar);
            isElemCohesiveFull       = any(ismember(edgesInElem,listEdgeCohesiveFull),2);
            edgesInCohElemReal = edgesInElem(isElemCohesiveFull,:);

            temp             = ismember(edgesInCohElemReal, listEdgeCohesiveFull);
            cohElemToEdge    = sum(edgesInCohElemReal .* temp, 2);    % (nElemCoh x 1) 
            centerElem       = centerElemsInCohesiveEdge;         % (nElemCoh x ndim)
            centerEdge       = centerEdges(cohElemToEdge,:);      % (nElemCoh x ndim)
            vectorEdgeToElem = centerElem - centerEdge;

            temp = zeros(size(obj.baseMesh.edges.nodesInEdges,1),2);


            temp(listEdgeCohesiveFull,:) = normalsFull;
            normalsFull = temp(cohElemToEdge,:);

            dotProduct  = sum(vectorEdgeToElem.*normalsFull,2);
            signs = sign(dotProduct);
            signs = signs > 0;
            isLeft = logical(signs);
            isRight = not(isLeft);
        end

        function newCoord = duplicateNodes(obj)
            newCoord    = obj.baseMesh.coord;
            listEdgeCohesiveFull = find(obj.isEdgeCohesiveAuxiliar | obj.isEdgeCohesiveReal);
            nodesInEdgesCohesiveFull = obj.baseMesh.edges.nodesInEdges(listEdgeCohesiveFull,:);
            listNodeCohesiveFull = unique(nodesInEdgesCohesiveFull');

            % nNodeCohesive    = length(obj.listNodeCohesive);
            isNodeCohesiveFull   = false(size(obj.baseMesh.coord,1),1);
            isNodeCohesiveFull(listNodeCohesiveFull) = true;
            nNodeCohesiveFull = sum(isNodeCohesiveFull);

            duplicated  = newCoord(isNodeCohesiveFull, :);
            newCoord    = [newCoord; duplicated];
            obj.pairsMatrixFull = [sort(listNodeCohesiveFull) , linspace(obj.baseMesh.nnodes+1, ...
                obj.baseMesh.nnodes+nNodeCohesiveFull, nNodeCohesiveFull)'];
            col1 = obj.baseMesh.edges.nodesInEdges(obj.isEdgeCohesiveReal,:)';
            col2 = obj.getPair(obj.pairsMatrixFull, col1);
            obj.pairsMatrixReal = [col1,col2];
        end

        function newConnec = updateConnecOfLeftElements(obj,isLeft,newCoord)
            listEdgeCohesiveReal = find(obj.isEdgeCohesiveReal);
            isElemCohesiveReal      = any(ismember(obj.baseMesh.edges.edgesInElem,listEdgeCohesiveReal),2);
            listElemNextCohesiveReal = find(isElemCohesiveReal);

            listEdgeCohesiveFull = find(obj.isEdgeCohesiveReal | obj.isEdgeCohesiveAuxiliar);
            isElemCohesiveFull      = any(ismember(obj.baseMesh.edges.edgesInElem,listEdgeCohesiveFull),2);
            listElemNextCohesiveFull = find(isElemCohesiveFull);

            listLeftElems  = listElemNextCohesiveFull(isLeft);

            % nodes = obj.pairsMatrixFull(:,1);
            % coords = obj.baseMesh.coord(nodes);

            cohesiveConnec = [obj.pairsMatrixReal(1:end-1,1), obj.pairsMatrixReal(2:end,1),...
                obj.pairsMatrixReal(2:end,2), obj.pairsMatrixReal(1:end-1,2)];
            connec         = [obj.baseMesh.connec; cohesiveConnec];

            % listCohesiveElems = ((size(connec,1)-size(cohesiveConnec,1)+1):size(connec,1))';
            oldLeftConnec = connec(listLeftElems,:);
            newLeftConnec      = oldLeftConnec;

            idx = ismember(oldLeftConnec,obj.pairsMatrixReal(:,1));
            newLeftConnec(idx) = arrayfun(@(x) obj.getPair(obj.pairsMatrixReal,x), oldLeftConnec(idx));

            connec(listLeftElems,:) = newLeftConnec;
            newConnec = connec;
        end

        function pair = getPair(~,pairMatrix,n)
            pair = pairMatrix(ismember(pairMatrix(:,1),n),2);
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
            
            rightNodes   = obj.pairsMatrixFull(:,1);
            leftNodes    = obj.getPair(obj.pairsMatrixFull,rightNodes);
            newCoord(rightNodes,:) = newCoord(rightNodes,:) - shiftVector*obj.separation/2;
            newCoord(leftNodes,:)  = newCoord(leftNodes,:) + shiftVector*obj.separation/2;   
        end

        function createLineMesh(obj)
            coord     = obj.fullMesh.coord;
            nMidNodes = size(obj.pairsMatrixReal,1);
            midCoord= (obj.pairsMatrixReal(:,1)+ ...
                coord(obj.pairsMatrixReal(:,2),:))./2;      
            midConnec =  [(1:nMidNodes-1)' (2:nMidNodes)'];
            s.connec = midConnec;
            s.coord  = midCoord;
            s.kFace  = -1;
            obj.mesh = Mesh.create(s);
        end


        function  newConnec = fixFinalElement(obj,cParams,newConnec)
            centerEdges      = obj.computeCenterEdge();
            listEdgeCohesiveFull = find(cParams.isFracturedLine(centerEdges));
            difference = listEdgeCohesiveFull(not(ismember(listEdgeCohesiveFull,obj.listEdgeCohesiveReal)));
            if not(isempty(difference))
                nodes = obj.baseMesh.edges.nodesInEdges(difference,:);
                nodes = unique(nodes);
                changedNode = obj.listNodeCohesive(ismember(obj.listNodeCohesive,nodes));
                idx = sum(ismember(obj.baseMesh.edges.nodesInEdges(difference,:),changedNode),2) == 1;
                uniqueEdge = difference(idx);
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
end
