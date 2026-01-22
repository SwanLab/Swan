classdef CohesiveMesh < handle
    
    properties (Access = public)
        baseMesh 
        mesh

        isNodeCohesive
        isElemCohesive
        isEdgeCohesive

        newCoord
        newConnec

        listNodeCohesive
        listElemCohesive
        listEdgeCohesive

        nNodeCohesive
        nElemCohesive
        pairsMatrix

        separation

        tangents
        normals
        
        cohElem2Edge


        elemsInCohesiveEdge
        centerElemsInCohesiveEdge
        centerEdges

        isLeft
        isRight


        % baseMesh i crear nova malla --> Mesh.create(s)
        % arreglar 1x1
        % cond contorn
        % triangle
        % logical

    end
    
    properties (Access = private)

    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveMesh()
            
            obj.init()
            obj.baseMeshCreator(2)
            
            obj.detectFracturedEdges()
                % fer lo del boundary amb els punts mitjos dels edges


            obj.computeCenterElements %(of fractureEdge)
            obj.computenormals()
            obj.computeIsLeftIsRight()
            obj.duplicator()
            obj.updateConnecOfLeftElements()
            % shiftCoordOfLeftAndRightElements();
            

            % obj.updateConnec()
            % obj.newMesh()
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

            obj.isElemCohesive = any(reshape(obj.isNodeCohesive(obj.baseMesh.connec), size(obj.baseMesh.connec)),2);
            obj.listElemCohesive = find(obj.isElemCohesive);

            edgesInElem = obj.baseMesh.edges.edgesInElem;
            edgesInCohElem = edgesInElem(obj.isElemCohesive,:);

            temp = ismember(edgesInCohElem, obj.listEdgeCohesive);
            obj.cohElem2Edge = sum(edgesInCohElem .* temp, 2);

        end

        function computeCenterElements(obj)

            bariCenters = obj.baseMesh.computeBaricenter';
            obj.centerElemsInCohesiveEdge = bariCenters(obj.listElemCohesive,:); %Ordenats segons obj.listElemCohesive

            % edgesInElem = obj.baseMesh.edges.edgesInElem;
            % 
            % nElem  = size(edgesInElem,1);
            % nEdges = max(edgesInElem(:));
            % 
            % E = repelem((1:nElem)', size(edgesInElem,2));  % NOT nEdgeByElem
            % temp = edgesInElem';
            % G = temp(:);
            % 
            % edgeElemMat = sparse(G, E, true, nEdges, nElem);
            % 
            % 
            % 
            % 
            % 
            %                 edgesInCohElem = edgesInElem(obj.isElemCohesive,:);
            % 
            %                 nElem  = size(edgesInCohElem,1);
            %                 nEdges = max(edgesInCohElem(:));
            % 
            %                 E = repelem((1:nElem)', size(edgesInCohElem,2));  % NOT nEdgeByElem
            %                 temp = edgesInCohElem';
            %                 G = temp(:);
            % 
            %                 edgeCohElemMat = sparse(G, E, true, nEdges, nElem);
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % 
            % for i=1:length(obj.listEdgeCohesive)
            %     obj.elemsInCohesiveEdge{i} = find(edgeElemMat(obj.listEdgeCohesive(i),:));
            %     % obj.centerElemsInCohesiveEdge{i} = num2cell(bariCenters(:, obj.elemsInCohesiveEdge{i}), 1);
            % end
            % % obj.centerElemsInCohesiveEdge=obj.centerElemsInCohesiveEdge';
            % 
            % %centerElemsInCohesiveEdge{cohesive edge}{element(1 o 2)} es una cell array, solament
            % %apareixen les coordenades dels elements cohesius, ordenats per
            % %edges cohesius, contempla la possibilitat de que un edge
            % %tingui mes de 1 elements
            % 
            % 
            % 



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

            obj.newConnec = obj.baseMesh.connec;
            cohesiveConnec = [obj.pairsMatrix(1:end-1,1), obj.pairsMatrix(2:end,1), obj.pairsMatrix(2:end,2), obj.pairsMatrix(1:end-1,2)];
            obj.newConnec = [obj.newConnec; cohesiveConnec];

        end

        function updateConnecOfLeftElements(obj)

            listLeftElems = obj.listElemCohesive(obj.isLeft);
            nodesInEdges = obj.baseMesh.edges.nodesInEdges; 

            for i =1:length(listLeftElems)
                e = listLeftElems(i);

                edge = obj.cohElem2Edge(e);
                replacedNodes = nodesInEdges(edge,:)';

                newNodes = obj.getPair(replacedNodes)';

                idx = ismember(obj.newConnec(e,:),replacedNodes');
                obj.newConnec(e,idx) = newNodes;
            end

        end

        function pair = getPair(obj,n)
            pair = obj.pairsMatrix(ismember(obj.pairsMatrix(:,1), n), 2 );

        end
    
        function newMesh(obj)
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
    
    
    end
        

end
