classdef PathToBoundaryComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       mesh 
       singularityCoord
       boundaryPointCoord       
    end
    
    properties (Access = private)
       pathVector
       closestVertex
       boundaryNodes
       connectedVertex       
       pathVertices       
    end
    
    methods (Access = public)
        
        function obj = PathToBoundaryComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.createStraightPathVector();
            obj.computeClosestVertex();
            obj.computeBoundaryNodes();
            obj.computeConnectedVertex();
            obj.computePathVertices();
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
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.singularityCoord   = cParams.singularityCoord;
            obj.boundaryPointCoord = cParams.boundaryPointCoord;
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
        
        function computeConnectedVertex(obj)
            obj.mesh.computeEdges();
            s.edges = obj.mesh.edges;
            c = ConnectedVertexesComputer(s);
            obj.connectedVertex = c;
        end
        
        function itIs = isInBoundary(obj,node)
            nodesB = obj.boundaryNodes;
            itIs = any(node == nodesB);
        end
        
        function computePathVertices(obj)            
            i = 1;
            vertex = obj.closestVertex;
            allVertices(i) = vertex;            
            while ~obj.isInBoundary(vertex)
                vertex = obj.selectNewOptimalVertex(vertex);
                i = i +1;  
                allVertices(i) = vertex;                                
            end
            obj.pathVertices = allVertices;
        end
        
        function newVertex = selectNewOptimalVertex(obj,oldVertex)
            tVertexes = obj.computeTrialVertexes(oldVertex);
            newVertex = obj.computeOptimalVertex(oldVertex,tVertexes);   
        end
        
        function tVertexes = computeTrialVertexes(obj,vertex)
            c = obj.connectedVertex;
            tVertexes = c.compute(vertex);
        end             
        
        function newVertexes = computeOptimalVertex(obj,currentVertex,trialVertexes)             
            s.currentVertexCoord  = obj.mesh.coord(currentVertex,:);
            s.trialVertexesCoord  = obj.mesh.coord(trialVertexes,:);
            s.boundaryVertexCoord = obj.boundaryPointCoord;
            s.lineVector          = obj.pathVector;
            n = NextVertexToBoundaryComputer(s);
            iD = n.compute();
            newVertexes = trialVertexes(iD);            
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
            cV = obj.pathVertices;
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            plot(x,y,'g-')
        end
        
    end
    
    
end