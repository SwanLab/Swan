classdef PathToBoundaryComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
       mesh
       coordSingularity
       coordFinal
       pathVector
       closestVertex
       vert
    end
    
    methods (Access = public)
        
        function obj = PathToBoundaryComputer()
            obj.createMesh();
            obj.createSingularityPoint();
            obj.createFinalPathPoint();
            obj.createStraightPathVector();
            obj.computeClosestVertex();
            obj.computePath();
            obj.plot();
        end

        function plot(obj)
            figure()
            hold on
            obj.plotMesh();
            obj.plotStraightPath();
            obj.plotClosestVertex();
            obj.plotVerticesPath();
        end        
        
    end
    
    methods (Access = private)
        
        function createMesh(obj)
            xmin = 0;
            xmax = 2;
            ymin = 0;
            ymax = 3;
            nP = 20000;           
            nB = 100;            
            xb = xmin + (xmax - xmin)*([0;rand(nB-2,1);1]);
            yb = xmin + (xmax - xmin)*([0;rand(nB-2,1);1]);
            

            xI = xmin + (xmax - xmin)*(rand(nP,1));
            yI = ymin + (ymax - ymin)*(rand(nP,1));            

            xT = [xb;xmin*ones(nB,1);xb;xmax*ones(nB,1);xI];
            yT = [ymin*ones(nB,1);yb;ymax*ones(nB,1);yb;yI];

            s.coord(:,1) = xT;
            s.coord(:,2) = yT;
            s.connec = delaunay(s.coord);
            
            xv = linspace(xmin,xmax,3);
            yv = linspace(ymin,ymax,4);
            [X,Y] = meshgrid(xv,yv);
            [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
             s.coord  = V(:,1:2);
             s.connec = F;    

            m = Mesh(s);
            obj.mesh = m;
        end

        function createSingularityPoint(obj)
            obj.coordSingularity = [0.59 2.7];
        end

        function createFinalPathPoint(obj)
            obj.coordFinal = [0.59 0];
        end

        function createStraightPathVector(obj)
            OA = obj.coordSingularity;            
            OB = obj.coordFinal;
            AB = OB - OA;
            obj.pathVector = AB ;
        end        

        function computeClosestVertex(obj)
            xy  = obj.mesh.coord;
            xyS = obj.coordSingularity;
            xD  = xy(:,1) - xyS(:,1);
            yD  = xy(:,2) - xyS(:,2);
            dist = sqrt(xD.^2 + yD.^2);
            [~,iD] = min(dist);
            obj.closestVertex = iD;
        end

        function computePath(obj)
            nodesB = boundary(obj.mesh.coord);
            obj.mesh.computeEdges();

            
            isInBoundary = false;
            i = 1;
            currentVertex = obj.closestVertex;
            vertices(i) = currentVertex;
            while ~isInBoundary
            nodeA = 1;
            nodeB = 2;            
            [candA,maxScA,pA] = obj.computeVertexV(currentVertex,nodeA,nodeB);


            nodeA = 2;
            nodeB = 1;            
            [candB,maxScB,pB] = obj.computeVertexV(currentVertex,nodeA,nodeB);

            tot = [maxScA, maxScB];
            pT = [pA,pB];
            idT = [candA,candB];
            [~,iDt] = obj.selectVertex(tot,pT);
            newCurrentVertex = idT(iDt);

            isInBoundary = any(newCurrentVertex == nodesB);
            currentVertex = newCurrentVertex;
            i = i +1;
            vertices(i) = newCurrentVertex;
            end
            obj.vert = vertices;
        end
        
        function [candB,csTmax,pB] = computeVertexV(obj,currentVertex,nodeA,nodeB)
            edges = obj.mesh.edges;
            edgeTrial   = find(edges.nodesInEdges(:,nodeA) == currentVertex);
            vertexB = edges.nodesInEdges(edgeTrial,nodeB);
            coordCurrentVertex = obj.mesh.coord(currentVertex,:);
            coordVertexB = obj.mesh.coord(vertexB,:);
            vectorTrial(:,1) = coordVertexB(:,1)  - coordCurrentVertex(1);
            vectorTrial(:,2) = coordVertexB(:,2)  - coordCurrentVertex(2);            
            direction = obj.pathVector;
            csT = obj.computeScalarProduct(vectorTrial,direction);
            %directionEnd(:,1)   = obj.coordFinal(:,1) - coordVertexB(:,1);             
            %directionEnd(:,2)   = obj.coordFinal(:,2) - coordVertexB(:,2); 
            directionEnd(:,1)   = obj.coordFinal(:,1) - coordCurrentVertex(:,1);                         
            directionEnd(:,2)   = obj.coordFinal(:,2) - coordCurrentVertex(:,2); 
            p = obj.computeProjection(obj.pathVector,vectorTrial,directionEnd);
            [maxScB,iD] = obj.selectVertex(csT,p);
            candB = vertexB(iD);
            pB = p(iD);
            csTmax = csT(iD);
        end

        function csT = computeScalarProduct(obj,a,b)
            scalarProdVD = obj.scalarProduct(a,b);
            scalarProdVV = obj.scalarProduct(a,a);
            scalarProdDD = obj.scalarProduct(b,b);
            csT = scalarProdVD./(sqrt(scalarProdVV).*sqrt(scalarProdDD));
        end

        function p = computeProjection(obj,AB,CE,CB)
            w = obj.computeUnitVector(AB);
            DB = obj.scalarProduct(CB,w)*w;
            CD = CB - DB;
            u = obj.computeUnitVector(CD);
            p = obj.scalarProduct(CE,u);
        end

        function u = computeUnitVector(obj,AB)
            sqAB = obj.scalarProduct(AB,AB);
            nAB = sqrt(sqAB);
            u(:,1) = AB(:,1)./nAB;
            u(:,2) = AB(:,2)./nAB;
        end

        function [maxScB,iD] = selectVertex(obj,csT,p)
            isP = csT > 0;
            p(~isP) = -inf;
            [maxScB,iD] = max(p);
            %[maxScB,iD] = max(csT);
        end

        function ab = scalarProduct(obj,a,b)
            ab = a(:,1).*b(:,1) + a(:,2).*b(:,2); 
        end

        function plotMesh(obj)
            obj.mesh.plot();
        end

        function plotStraightPath(obj)
            xP = [obj.coordSingularity(:,1),obj.coordFinal(:,1)];
            yP = [obj.coordSingularity(:,2),obj.coordFinal(:,2)];
            plot(xP,yP,'+-')            
        end

        function plotClosestVertex(obj)
            cV = obj.closestVertex;            
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            plot(x,y,'r+')
        end

        function plotVerticesPath(obj)
            cV = obj.vert;            
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            plot(x,y,'g-')
        end


    end
    
end