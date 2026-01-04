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

        nnodeCohesive
        nElemCohesive
        pairsMatrix

        separation


        
        
        centerEdges

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
            obj.baseMeshCreator(3)
            
            obj.detectFracturedEdges()
                % fer lo del boundary amb els punts mitjos dels edges


            obj.computeCenterElements %(of fractureEdge)
            % computeNormals()
            % computeIsLeftIsRight()
            % addDuplicateCoords();
            % updateConnecOfLeftElements()
            % shiftCoordOfLeftAndRightElements();
            
            % obj.duplicator()
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
            ymin = min(obj.baseMesh.coord(:,2));
            obj.centerEdges=obj.computeCenterEdge;
            obj.isEdgeCohesive = abs(obj.centerEdges(:,2)) == ymin;
            obj.listEdgeCohesive = find(obj.isEdgeCohesive==1);

            obj.isNodeCohesive = abs(obj.baseMesh.coord(:,2)) == ymin;
            obj.listNodeCohesive = find(obj.isNodeCohesive == 1);

            obj.isElemCohesive = any(reshape(obj.isNodeCohesive(obj.baseMesh.connec), size(obj.baseMesh.connec)),2);
        end

        function computeCenterElements(obj)
            edgesInElem = obj.baseMesh.edges.edgesInElem;
            nEdges = max(edgesInElem(:));
               
            E = repelem((1:obj.baseMesh.nelem)', obj.baseMesh.edges.nEdgeByElem);   
            G = edgesInElem(:);

            edgeElemMat = sparse(G, E, true, nEdges, obj.baseMesh.nelem);
            
            for e = obj.listEdgeCohesive
                elemsInCohesiveEdge{e} = find(edgeElemMat(e,:));
            end

            for i=1:length(obj.listEdgeCohesive)
                centerElemsInCohesiveEdge{i}=elemsInCohesiveEdge{i}(:)

            end
        end



        function duplicator(obj)

            % ymin = min(obj.baseMesh.coord(:,2));
            % obj.isCohesive = abs(obj.baseMesh.coord(:,2)) == ymin; 
            % obj.nnodeCohesive = sum(obj.isCohesive) * 2;
            % obj.listNodeCohesive = find(obj.isCohesive == 1);

            obj.newCoord = obj.baseMesh.coord;

            dispV      = [0, 1];
            duplicated = obj.newCoord(obj.isNodeCohesive, :) + obj.separation*dispV;
            obj.newCoord = [obj.newCoord; duplicated];

            obj.pairsMatrix = [obj.listNodeCohesive , (size(obj.isNodeCohesive,1)+1:1:obj.nnodeCohesive/2+size(obj.isNodeCohesive,1))'];
                obj.pairsMatrix = [obj.pairsMatrix , dispV(1)*ones(obj.nnodeCohesive/2,1), dispV(2)*ones(obj.nnodeCohesive/2,1)];

            obj.isNodeCohesive = [obj.isNodeCohesive; ones(obj.nnodeCohesive/2,1)];

                % scatter(obj.newCoord(:,1), obj.newCoord(:,2), 'filled');

        end

        function updateConnec(obj)
            obj.newConnec = obj.baseMesh.connec;

            isCohesiveElem = any(reshape(obj.isNodeCohesive(obj.newConnec), size(obj.newConnec)),2);
            obj.listElemCohesive = find(isCohesiveElem);

            % for i = obj.listNodeCohesive'
            %     for n = 1:4
            %         if obj.isCohesive(obj.newConnec(i,n))
            %             if n == 4
            %                 vec1 = obj.newCoord(obj.newConnec(i,4),:) - obj.newCoord(obj.newConnec(i,1),:);
            %             else
            %                 vec1 = obj.newCoord(obj.newConnec(i,n),:) - obj.newCoord(obj.newConnec(i,n+1),:);
            %             end
            %         end
            %         idx = find(obj.pairsMatrix(:,1) == obj.newConnec(i,n));
            %         vec2 = obj.pairsMatrix(idx,3:end);
            %         if dot(vec1',vec2') >0
            %             obj.newConnec(i,n) = obj.pairsMatrix(idx,2);
            %         end
            %     end
            % end

            for e = obj.listElemCohesive'
                nodes = obj.newConnec(e,:)';
                coords = obj.newCoord(nodes,:);
                
                cohesiveNodes = nodes((obj.isNodeCohesive(nodes)==1));
                originalArea = polyarea(coords(:,1), coords(:,2));
                
                
                for cn = cohesiveNodes'
                    pair = getPair(obj,cn);

                    checkNodes = nodes;
                    checkNodes(checkNodes == cn) = pair;
                    checkCoords = obj.newCoord(checkNodes,:);

                    newArea = polyarea(checkCoords(:,1), checkCoords(:,2));

                    if newArea<originalArea
                        obj.newConnec(e,find(obj.newConnec(e,:) == cn)) = pair;
                    end

                end

            end

            cohesiveConnec = [obj.pairsMatrix(1:end-1,1), obj.pairsMatrix(2:end,1), obj.pairsMatrix(2:end,2), obj.pairsMatrix(1:end-1,2)];
            obj.newConnec = [obj.newConnec; cohesiveConnec];
            obj.nElemCohesive = size(cohesiveConnec,1);
            
            obj.baseMesh.coord = obj.newCoord;
            obj.baseMesh.connec = obj.newConnec;

        end

        function pair = getPair(obj,n)
            pair = obj.pairsMatrix(obj.pairsMatrix(:,1) == n, 2);
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
    
