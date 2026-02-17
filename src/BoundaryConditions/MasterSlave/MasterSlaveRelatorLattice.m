classdef MasterSlaveRelatorLattice < handle
    
    properties
        coord
        a1
        a2
        tol = 1e-8
        msRelation
        meshType  
    end
    
    methods
        function obj = MasterSlaveRelatorLattice(coord, a1, a2, varargin)
            obj.coord = coord;
            obj.a1 = a1(:);
            obj.a2 = a2(:);
            
            
           
            obj.detectMeshType();
            
            
            for i = 1:2:length(varargin)
                switch varargin{i}
                    case 'meshType'
                        obj.meshType = varargin{i+1};
                end
            end
            
            obj.computeRelation();
        end
        
        function detectMeshType(obj)
            
            angle = acos(dot(obj.a1, obj.a2)/(norm(obj.a1)*norm(obj.a2))) * 180/pi;
            
            if abs(angle - 90) < 1
                if abs(norm(obj.a1) - norm(obj.a2)) < 1e-6
                    obj.meshType = 'Square';
                else
                    obj.meshType = 'Rectangle';
                end
            elseif abs(angle - 60) < 1 || abs(angle - 120) < 1
                obj.meshType = 'Hexagon';
            else
                obj.meshType = 'General';
            end
        end
        
        function computeRelation(obj)
            switch obj.meshType
                case {'Square', 'Rectangle', 'General'}
                    
                    obj.computeRelationParallelogram();
                case 'Hexagon'
                    
                    obj.computeRelationHexagon();
                otherwise
                    obj.computeRelationParallelogram();
            end
        end
        
        function computeRelationParallelogram(obj)
            
            A = [obj.a1 obj.a2];
            invA = inv(A);
            
            local = (invA * obj.coord')';
            xi  = local(:,1);
            eta = local(:,2);
            
            
            isXiMin = abs(xi - min(xi)) < obj.tol;
            isXiMax = abs(xi - max(xi)) < obj.tol;
            isEtaMin = abs(eta - min(eta)) < obj.tol;
            isEtaMax = abs(eta - max(eta)) < obj.tol;
            
            
            corners = find((isXiMin & isEtaMin) | (isXiMin & isEtaMax) | ...
                          (isXiMax & isEtaMin) | (isXiMax & isEtaMax));
            isXiMin(corners) = false; isXiMax(corners) = false;
            isEtaMin(corners) = false; isEtaMax(corners) = false;
            
            
            pairs = [];
            pairs = [pairs; obj.matchFaces(isXiMin, isXiMax, eta, true)];
            pairs = [pairs; obj.matchFaces(isEtaMin, isEtaMax, xi, true)];
            
            obj.msRelation = pairs;
        end
        
        function computeRelationHexagon(obj)
            
            vertices = obj.findHexagonVertices();
            nVertices = length(vertices);
            
            
            vertices = obj.orderVerticesCCW(vertices);
            
            
            pairs = [];
            for i = 1:3
                v1 = vertices(i);
                v2 = vertices(mod(i, nVertices) + 1);
                v3 = vertices(i+3);
                v4 = vertices(mod(i+3, nVertices) + 1);
                
                
                nodesSide1 = obj.findNodesOnSide(v1, v2);
                nodesSide2 = obj.findNodesOnSide(v3, v4);
                
                
                nodesSide1 = setdiff(nodesSide1, [v1, v2]);
                nodesSide2 = setdiff(nodesSide2, [v3, v4]);
                
               
                nodesSide1 = obj.sortNodesAlongSide(nodesSide1, v1, v2);
                nodesSide2 = obj.sortNodesAlongSide(nodesSide2, v3, v4);
                
                
                nodesSide2 = flipud(nodesSide2);
                
                
                n = min(length(nodesSide1), length(nodesSide2));
                pairs = [pairs; nodesSide1(1:n), nodesSide2(1:n)];
            end
            
            obj.msRelation = pairs;
        end
        
        function pairs = matchFaces(obj, faceA, faceB, coordSort, invertSlave)
            
            nodesA = find(faceA);
            nodesB = find(faceB);
            
            [~, iA] = sort(coordSort(nodesA));
            [~, iB] = sort(coordSort(nodesB));
            
            nodesA = nodesA(iA);
            nodesB = nodesB(iB);
            
            
            if invertSlave
                nodesB = flipud(nodesB);
            end
            
            n = min(length(nodesA), length(nodesB));
            pairs = [nodesA(1:n), nodesB(1:n)];
        end
        
        
        function vertices = findHexagonVertices(obj)
            
            v{1} = [0, 0];
            v{2} = obj.a1';
            v{3} = obj.a2';
            v{4} = obj.a1' + obj.a2';
            v{5} = 2*obj.a1' + obj.a2';
            v{6} = obj.a1' + 2*obj.a2';
            
            
            vertices = zeros(6,1);
            for i = 1:6
                dist = sqrt(sum((obj.coord - v{i}).^2, 2));
                [~, vertices(i)] = min(dist);
            end
        end
        
        function vertices = orderVerticesCCW(obj, vertices)
            
            centroid = mean(obj.coord(vertices, :), 1);
            vectors = obj.coord(vertices, :) - centroid;
            angles = atan2(vectors(:,2), vectors(:,1));
            [~, idx] = sort(angles);
            vertices = vertices(idx);
        end
        
        function nodes = findNodesOnSide(obj, v1, v2)
            
            lineVec = obj.coord(v2,:) - obj.coord(v1,:);
            lineNorm = lineVec / norm(lineVec);
            
            
            vecToNode = obj.coord - obj.coord(v1,:);
            perpDist = abs(vecToNode(:,1)*lineNorm(2) - vecToNode(:,2)*lineNorm(1));
            
            
            onLine = perpDist < obj.tol;
            
            
            proj = vecToNode * lineNorm';
            between = proj >= -obj.tol & proj <= norm(lineVec) + obj.tol;
            
            nodes = find(onLine & between);
        end
        
        function sortedNodes = sortNodesAlongSide(obj, nodes, v1, v2)
            
            vec = obj.coord(v2,:) - obj.coord(v1,:);
            proj = (obj.coord(nodes,:) - obj.coord(v1,:)) * (vec/norm(vec))';
            [~, idx] = sort(proj);
            sortedNodes = nodes(idx);
        end
        
        function r = getRelation(obj)
            r = obj.msRelation;
        end
        
        function plotPairs(obj)
           
            figure; hold on;
            plot(obj.coord(:,1), obj.coord(:,2), 'b.');
            
            for i = 1:size(obj.msRelation, 1)
                master = obj.msRelation(i, 1);
                slave = obj.msRelation(i, 2);
                
                plot(obj.coord(master,1), obj.coord(master,2), 'go', 'MarkerSize', 8);
                plot(obj.coord(slave,1), obj.coord(slave,2), 'ro', 'MarkerSize', 8);
                
                
                plot([obj.coord(master,1), obj.coord(slave,1)], ...
                     [obj.coord(master,2), obj.coord(slave,2)], 'k-');
            end
            axis equal;
            title(sprintf('Master-Slave Pairs (%s mesh)', obj.meshType));
        end
    end
end