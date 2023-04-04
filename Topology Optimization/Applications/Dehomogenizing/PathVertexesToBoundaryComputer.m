classdef PathVertexesToBoundaryComputer < handle
    
    properties (Access = private)
       mesh 
       singularityCoord
    end
    
    properties (Access = private)
       boundaryPointCoord               
       pathVector
       closestVertex
       boundaryNodes
       connectedVertex     
       pathVertexes
       pathCells
    end
    
    methods (Access = public)
        
        function obj = PathVertexesToBoundaryComputer(cParams)
            obj.init(cParams)
        end
        
        function v = compute(obj)
            obj.computeBoundaryPointCoord();
            obj.createStraightPathVector();
            obj.computeClosestVertex();
            obj.computeBoundaryNodes();
            obj.computeMeshEdges();
            obj.computePath();
            v = obj.pathVertexes;          
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
              
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.singularityCoord   = cParams.singularityCoord;
        end
        
        function computeBoundaryPointCoord(obj)
        %    obj.computeToyBoundaryPoint();
         %  obj.computeBenchmarkBoundaryPoint();
            sC = obj.singularityCoord;
            bM = obj.mesh.createBoundaryMesh();
            for iB = 1:numel(bM)
               bmI = bM{iB};
               bC  = bmI.mesh.coord;
               dis = (bC(:,1) - sC(1)).^2 + (bC(:,2) - sC(2)).^2;
               [dM,im] = min(dis);
               nodes(iB) = bmI.globalConnec(im);
               dist(iB) = dM;
            end
            [~,iM] = min(dist);
            node = nodes(iM);
            cP = obj.mesh.coord(node,:);
           obj.boundaryPointCoord = cP;
            %obj.boundaryPointCoord(:,2)
        end
        
        function computeToyBoundaryPoint(obj)
            obj.boundaryPointCoord = [0 150];            
        end        

        function computeBenchmarkBoundaryPoint(obj)
           obj.boundaryPointCoord(:,2) = obj.singularityCoord(:,2);                       
           obj.boundaryPointCoord(:,1) = max(obj.mesh.coord(:,1));            
        end               

        function createStraightPathVector(obj)
            OA = obj.singularityCoord;            
            OB = obj.boundaryPointCoord;
            AB = OB - OA;
            obj.pathVector = AB ;
        end        

        function computeClosestVertex(obj)
            xy  = obj.mesh.coord;
            xyS = obj.singularityCoord;
            xD  = xy(:,1) - xyS(:,1);
            yD  = xy(:,2) - xyS(:,2);
            dist = sqrt(xD.^2 + yD.^2);
            [~,iD] = min(dist);
            obj.closestVertex = iD;
        end
        
        function computeBoundaryNodes(obj)
            coord  = obj.mesh.coord;
            nodesB = boundary(coord);            
            obj.boundaryNodes = nodesB;
        end
        
        function computeMeshEdges(obj)
            obj.mesh.computeEdges();            
        end

        function itIs = isInBoundary(obj,node)
            nodesB = obj.boundaryNodes;
            itIs = any(node == nodesB);
        end
        
        function computePath(obj)            
            i = 1;
            vertex       = obj.closestVertex;
            pVertexes    = vertex;
            while ~obj.isInBoundary(vertex)
                otherVertex = obj.mesh.computeConnectedVertex(vertex);
                iD          = obj.computeOptimalVertex(vertex,otherVertex);  
                newVertex   = otherVertex(iD);    
                i = i + 1;  
                vertex       = newVertex;                
                pVertexes(i) = vertex;              
            end
            obj.pathVertexes = pVertexes;
        end        
        
        function iD = computeOptimalVertex(obj,currentVertex,trialVertexes)             
            s.currentVertexCoord  = obj.mesh.coord(currentVertex,:);
            s.trialVertexesCoord  = obj.mesh.coord(trialVertexes,:);
            s.boundaryVertexCoord = obj.boundaryPointCoord;
            s.lineVector          = obj.pathVector;
            n = NextVertexToBoundaryComputer(s);
            iD = n.compute();
        end            
        
        function plotMesh(obj)
            obj.mesh.plot();
        end
        
        function plotStraightPath(obj)
            xP = [obj.singularityCoord(:,1),obj.boundaryPointCoord(:,1)];
            yP = [obj.singularityCoord(:,2),obj.boundaryPointCoord(:,2)];
            plot(xP,yP,'+-')
        end
        
        function plotClosestVertex(obj)
            cV = obj.closestVertex;
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            plot(x,y,'r+')
        end
        
        function plotVerticesPath(obj)
            cV = obj.pathVertexes;
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            plot(x,y,'g-','LineWidth',5)
        end
        
    end
    
    
    
end