classdef CohesiveMesh < handle
    
    properties (Access = public)
        baseMesh 
        mesh
        subMesh

        isNodeCohesive
        isElemCohesive
        isEdgeCohesive

        newCoord
        newConnec

        listNodeCohesive
        listElemNextCohesive
        listCohesiveElems
        listEdgeCohesive

        nNodeCohesive
        pairsMatrix

        separation

        tangents
        normals
        
        cohElem2Edge


        centerElemsInCohesiveEdge
        centerEdges

        isLeft
        isRight


    end
    
    properties (Access = private)

    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveMesh()
            
            obj.init()
            obj.baseMeshCreator(3)
            
            obj.detectFracturedEdges()
            obj.computeCenterElements 
            obj.computenormals()
            obj.computeIsLeftIsRight()
            obj.duplicator()
            obj.updateConnecOfLeftElements()
            obj.shiftCoordOfLeftAndRightElements();
            
            obj.newMesh()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.separation = 0.1;
        end
        
        function baseMeshCreator(obj,n)
            obj.baseMesh = UnitQuadMesh(n,n);
        end

        function detectFracturedEdges(obj)
            obj.centerEdges=obj.computeCenterEdge;

            ymin = min(obj.baseMesh.coord(:,2));
            obj.isEdgeCohesive = abs(obj.centerEdges(:,2)) == ymin;
            obj.listEdgeCohesive = find(obj.isEdgeCohesive);

            obj.isNodeCohesive = abs(obj.baseMesh.coord(:,2)) == ymin;
            obj.listNodeCohesive = find(obj.isNodeCohesive);
            obj.nNodeCohesive = sum(obj.isNodeCohesive);

            obj.isElemCohesive = any(ismember(obj.baseMesh.edges.edgesInElem,obj.listEdgeCohesive),2);
            obj.listElemNextCohesive = find(obj.isElemCohesive);

            edgesInElem = obj.baseMesh.edges.edgesInElem;
            edgesInCohElem = edgesInElem(obj.isElemCohesive,:);

            temp = ismember(edgesInCohElem, obj.listEdgeCohesive);
            obj.cohElem2Edge = sum(edgesInCohElem .* temp, 2);

        end

        function computeCenterElements(obj)
            
            bariCenters = obj.baseMesh.computeBaricenter';
            obj.centerElemsInCohesiveEdge = bariCenters(obj.listElemNextCohesive,:); %Ordenats segons obj.listElemCohesive

        end

        
        function computenormals(obj)
            nodesInEdges = obj.baseMesh.edges.nodesInEdges;

            coord = obj.baseMesh.coord;
            
            nodes = nodesInEdges(obj.listEdgeCohesive,:);
            
            coords1 = coord(nodes(:,1),:);   % nCohEdges x 2
            coords2 = coord(nodes(:,2),:);   % nCohEdges x 2
            
            obj.tangents = coords2 - coords1;        % nCohEdges x 2
            
            obj.normals = [obj.tangents(:,2),obj.tangents(:,1)];
        end


        function computeIsLeftIsRight(obj)
            edges = obj.cohElem2Edge;                 % (nElemCoh x 1)
            centerElem = obj.centerElemsInCohesiveEdge; % (nElemCoh x ndim)
            centerEdge = obj.centerEdges(edges,:);      % (nElemCoh x ndim)
            
            vectorEdgeToElem = centerElem - centerEdge;

            dotProduct = sum(vectorEdgeToElem.*obj.normals,2);
            signs = sign(dotProduct);
            obj.isLeft = logical(signs);
            obj.isRight = not(obj.isLeft);
        end


        function duplicator(obj)

            obj.newCoord = obj.baseMesh.coord;
            duplicated = obj.newCoord(obj.isNodeCohesive, :);
            obj.newCoord = [obj.newCoord; duplicated];

            obj.pairsMatrix = [obj.listNodeCohesive , linspace(obj.baseMesh.nnodes+1, ...
                obj.baseMesh.nnodes+obj.nNodeCohesive,obj.nNodeCohesive)'];


            obj.isNodeCohesive = [obj.isNodeCohesive; ones(obj.nNodeCohesive,1)];

        end

        function updateConnecOfLeftElements(obj)

            listLeftElems = obj.listElemNextCohesive(obj.isLeft);

            connec = obj.baseMesh.connec;
            cohesiveConnec = [obj.pairsMatrix(1:end-1,1), obj.pairsMatrix(2:end,1), obj.pairsMatrix(2:end,2), obj.pairsMatrix(1:end-1,2)];
            connec = [connec; cohesiveConnec];

            obj.listCohesiveElems = ((size(connec,1)-size(cohesiveConnec,1)+1):size(connec,1))';
            
            oldLeftConnec = connec(listLeftElems,:);
            idx = ismember(oldLeftConnec,obj.pairsMatrix(:,1));

            newLeftConnec = oldLeftConnec;
            newLeftConnec(idx) = arrayfun(@(x) obj.getPair(x), oldLeftConnec(idx));

            connec(listLeftElems,:) = newLeftConnec;
            
            
            obj.newConnec = connec;


        end

        function pair = getPair(obj,n)
            pair = obj.pairsMatrix(ismember(obj.pairsMatrix(:,1), n), 2 );

        end
    
        function newMesh(obj)
            s.connec = obj.newConnec;
            s.coord  = obj.newCoord;
            obj.mesh = Mesh.create(s);

            subplot(1,2,1)
            obj.baseMesh.plot;
            title('BaseMesh')
            
            subplot(1,2,2)
            obj.mesh.plot;
            title('CohesiveMesh')

        end
    
        function centerEdges = computeCenterEdge(obj)

            obj.baseMesh.computeEdges;
            nodes = obj.baseMesh.edges.nodesInEdges;    
            coord1 = obj.baseMesh.coord(nodes(:,1),:);
            coord2 = obj.baseMesh.coord(nodes(:,2),:);
            
            centerEdges = 0.5*(coord1 + coord2);

        end
        
        function shiftCoordOfLeftAndRightElements(obj)
                        
            normal = obj.normals;

            shiftVector = 0.5*([normal;zeros(1,size(normal,2))]+[zeros(1,size(normal,2));normal]);
            shiftVector(1,:) = normal(1,:); shiftVector(end,:) = normal(end,:);
            shiftVector = shiftVector./vecnorm(shiftVector,2,2); %esta ordenat segons listCohesiveNodes
            
            rightNodes = obj.listNodeCohesive;
            leftNodes = getPair(obj,rightNodes);

            obj.newCoord(rightNodes,:) = obj.newCoord(rightNodes,:) - shiftVector*obj.separation/2;
            obj.newCoord(leftNodes,:) = obj.newCoord(leftNodes,:) + shiftVector*obj.separation/2;
            
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
