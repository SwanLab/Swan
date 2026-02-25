classdef CohesiveMesh < handle
    
    properties (Access = public)
        baseMesh 
        mesh
        subMesh

        isNodeCohesive
        isElemCohesive

        listNodeCohesive
        listElemNextCohesive
        listCohesiveElems

        nNodeCohesive
        pairsMatrix

        edgesInCohElem

        separation

    end
    
    properties (Access = private)

    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveMesh(cParams)
            
            obj.init(cParams)
            obj.baseMeshCreator(3)
            
            listEdgeCohesive          = obj.detectFracturedEdges();
            centerElemsInCohesiveEdge = obj.computeCenterElements();
            normals                   = obj.computeNormals(listEdgeCohesive);
            [isLeft, isRight]         = obj.computeIsLeftIsRight(centerElemsInCohesiveEdge,normals);

            newCoord    = obj.duplicator();
            newConnec   = obj.updateConnecOfLeftElements(isLeft, isRight);
            newCoord    = obj.shiftCoordOfLeftAndRightElements(newCoord);
            
            obj.newMesh(newConnec, newCoord);
            obj.createSubMesh();
        end


        function plot(obj)
            
            subplot(1,2,1)
            obj.baseMesh.plot;
            title('BaseMesh')

            subplot(1,2,2)
            obj.mesh.plot;
            title('CohesiveMesh')
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.separation = 0.001;
            obj.baseMesh = cParams.baseMesh;
        end
        
        function baseMeshCreator(obj,n)
            obj.baseMesh = UnitQuadMesh(n,n);
        end

        function listEdgeCohesive = detectFracturedEdges(obj)
            
            centerEdges = computeCenterEdge();

            ymin = min(obj.baseMesh.coord(:,2));
            isEdgeCohesive = abs(centerEdges(:,2)) == ymin;
            listEdgeCohesive = find(isEdgeCohesive);

            obj.isNodeCohesive = abs(obj.baseMesh.coord(:,2)) == ymin;
            obj.listNodeCohesive = find(obj.isNodeCohesive);
            obj.nNodeCohesive = sum(obj.isNodeCohesive);

            obj.isElemCohesive = any(ismember(obj.baseMesh.edges.edgesInElem,listEdgeCohesive),2);
            obj.listElemNextCohesive = find(obj.isElemCohesive);

            edgesInElem = obj.baseMesh.edges.edgesInElem;
            obj.edgesInCohElem = edgesInElem(obj.isElemCohesive,:);



        end

        function centerElemsInCohesiveEdge = computeCenterElements(obj)
            
            bariCenters = obj.baseMesh.computeBaricenter';
            centerElemsInCohesiveEdge = bariCenters(obj.listElemNextCohesive,:); %Ordenats segons obj.listElemCohesive

        end

        
        function n = computeNormals(obj,listEdgeCohesive)
            nodesInEdges = obj.baseMesh.edges.nodesInEdges;

            coord = obj.baseMesh.coord;
            
            nodes = nodesInEdges(listEdgeCohesive,:);
            
            coords1 = coord(nodes(:,1),:);   % nCohEdges x 2
            coords2 = coord(nodes(:,2),:);   % nCohEdges x 2
            
            t = coords2 - coords1;        % nCohEdges x 2
            
            n = [t(:,2),t(:,1)];
        end


        function [isLeft, isRight] = computeIsLeftIsRight(obj,centerElemsInCohesiveEdge,normals)
            
            centerEdges=obj.computeCenterEdge;

            temp = ismember(obj.edgesInCohElem, listEdgeCohesive);
            cohElem2Edge = sum(obj.edgesInCohElem .* temp, 2); % (nElemCoh x 1)
                         
            centerElem = centerElemsInCohesiveEdge; % (nElemCoh x ndim)
            centerEdge = centerEdges(cohElem2Edge,:);      % (nElemCoh x ndim)
            
            vectorEdgeToElem = centerElem - centerEdge;

            dotProduct = sum(vectorEdgeToElem.*normals,2);
            signs = sign(dotProduct);
            isLeft = logical(signs);
            isRight = not(isLeft);
        end


        function newCoord = duplicator(obj)

            newCoord = obj.baseMesh.coord;
            duplicated = newCoord(obj.isNodeCohesive, :);
            newCoord = [newCoord; duplicated];

            obj.pairsMatrix = [obj.listNodeCohesive , linspace(obj.baseMesh.nnodes+1, ...
                obj.baseMesh.nnodes+obj.nNodeCohesive,obj.nNodeCohesive)'];


            obj.isNodeCohesive = [obj.isNodeCohesive; ones(obj.nNodeCohesive,1)];

        end

        function newConnec = updateConnecOfLeftElements(obj)

            listLeftElems = obj.listElemNextCohesive(isLeft);

            connec = obj.baseMesh.connec;
            cohesiveConnec = [obj.pairsMatrix(1:end-1,1), obj.pairsMatrix(2:end,1), obj.pairsMatrix(2:end,2), obj.pairsMatrix(1:end-1,2)];
            connec = [connec; cohesiveConnec];

            obj.listCohesiveElems = ((size(connec,1)-size(cohesiveConnec,1)+1):size(connec,1))';
            
            oldLeftConnec = connec(listLeftElems,:);
            idx = ismember(oldLeftConnec,obj.pairsMatrix(:,1));

            newLeftConnec = oldLeftConnec;
            newLeftConnec(idx) = arrayfun(@(x) obj.getPair(x), oldLeftConnec(idx));

            connec(listLeftElems,:) = newLeftConnec;
            newConnec = connec;

        end

        function pair = getPair(obj,n)
            pair = obj.pairsMatrix(ismember(obj.pairsMatrix(:,1), n), 2 );
        end
    
        function newMesh(obj,ne)
            s.connec = obj.newConnec;
            s.coord  = obj.newCoord;
            obj.mesh = Mesh.create(s);
        end
    
        function centerEdges = computeCenterEdge(obj)

            obj.baseMesh.computeEdges;
            nodes = obj.baseMesh.edges.nodesInEdges;    
            coord1 = obj.baseMesh.coord(nodes(:,1),:);
            coord2 = obj.baseMesh.coord(nodes(:,2),:);
            
            centerEdges = 0.5*(coord1 + coord2);

        end
        
        function newCoord = shiftCoordOfLeftAndRightElements(obj,newCoord,normals)

            shiftVector = 0.5*([normals;zeros(1,size(normals,2))]+[zeros(1,size(normals,2));normals]);
            shiftVector(1,:) = normals(1,:); shiftVector(end,:) = normals(end,:);
            shiftVector = shiftVector./vecnorm(shiftVector,2,2); %esta ordenat segons listCohesiveNodes
            
            rightNodes = obj.listNodeCohesive;
            leftNodes = getPair(obj,rightNodes);

            newCoord(rightNodes,:) = newCoord(rightNodes,:) - shiftVector*obj.separation/2;
            newCoord(leftNodes,:) = newCoord(leftNodes,:) + shiftVector*obj.separation/2;
            
        end


        function createSubMesh(obj)

            coord = obj.mesh.coord;
            nMidNodes = size(obj.pairsMatrix,1);
            
            midCoord= (coord(obj.listNodeCohesive,:)+ ...
                coord(obj.pairsMatrix(:,2),:))/2;
                        
            midConnec =  [(1:nMidNodes-1)' (2:nMidNodes)'];
            
                s.connec = midConnec;
                s.coord = midCoord;
                s.kFace = -1;
            obj.subMesh = Mesh.create(s);

        end

    
    end
        

end
