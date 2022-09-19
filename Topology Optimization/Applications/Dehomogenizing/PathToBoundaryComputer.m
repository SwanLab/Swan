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
       pathEdges
    end
    
    methods (Access = public)
        
        function obj = PathToBoundaryComputer(cParams)
            obj.init(cParams)
        end
        
        function [p,e] = compute(obj)
            obj.createStraightPathVector();
            obj.computeClosestVertex();
            obj.computeBoundaryNodes();
            obj.computeMeshEdges();
            obj.computePathVertices();
            p = obj.pathVertices;
            e = obj.pathEdges;
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
        
        function computeMeshEdges(obj)
            obj.mesh.computeEdges();            
        end

        function itIs = isInBoundary(obj,node)
            nodesB = obj.boundaryNodes;
            itIs = any(node == nodesB);
        end
        
        function computePathVertices(obj)            
            i = 1;
            vertex       = obj.closestVertex;
            pVertexes    = vertex;
            pEdges       = zeros(1);
            while ~obj.isInBoundary(vertex)
                edges       = obj.computeAllEdgesOfVertex(vertex);
                otherVertex = obj.computeOtherVertexOfEdge(edges,vertex);
                iD          = obj.computeOptimalVertex(vertex,otherVertex);  
                newVertex   = otherVertex(iD);    
                newEdge     = edges(iD);    
                pEdges(i)   = newEdge;
                i = i + 1;  
                vertex       = newVertex;                
                pVertexes(i) = vertex;
            end
            obj.pathVertices = pVertexes;
            obj.pathEdges    = pEdges;
        end
        

        function edges = computeAllEdgesOfVertex(obj,vertex)
            vertexInEdges = obj.mesh.edges.nodesInEdges;            
            isInEdge = vertexInEdges == vertex;            
            allEdges  = 1:size(isInEdge,1);
            allEdgesA(:,1) = allEdges(isInEdge(:,2));
            allEdgesB(:,1) = allEdges(isInEdge(:,1));
            edges = [allEdgesA;allEdgesB];            
        end
        
        function oV = computeOtherVertexOfEdge(obj,edge,vertex)
            vertexesInEdges = obj.mesh.edges.nodesInEdges;
            vertexesOfEdge  = vertexesInEdges(edge,:);
            oV = setdiff(vertexesOfEdge(:),vertex);
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
            cV = obj.pathVertices;
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            plot(x,y,'g-')
        end
        
    end
    
    
end