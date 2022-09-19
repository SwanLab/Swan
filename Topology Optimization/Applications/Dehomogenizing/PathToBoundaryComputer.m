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
       pathCells
    end
    
    methods (Access = public)
        
        function obj = PathToBoundaryComputer(cParams)
            obj.init(cParams)
        end
        
        function [v,e,c] = compute(obj)
            obj.createStraightPathVector();
            obj.computeClosestVertex();
            obj.computeBoundaryNodes();
            obj.computeMeshEdges();
            obj.computePath();
            v = obj.pathVertices;
            e = obj.pathEdges;
            c = obj.pathCells;
            
            obj.computeLeftRightElements();
            
            s.coord = obj.mesh.coord;
            s.connec = obj.mesh.connec(c,:);
            m = Mesh(s);
            m.plot()
            obj.plot();            
        end

        function plot(obj)
         %   figure()
            hold on
            obj.plotMesh();
            obj.plotStraightPath();
            obj.plotClosestVertex();
            obj.plotVerticesPath();
        end        
        
    end
    
    methods (Access = private)
        
        function computeLeftRightElements(obj)
            v = obj.pathVertices;
            %c = obj.pathCells;
            cellBar = obj.mesh.computeBaricenter()';
            cellRight = [];
            cellLeft = [];
            
            for i=2:(length(v)-1)
                vertexOld = v(i-1);
                vertexI   = v(i);
                vertexNew = v(i+1);
                coordOld = obj.mesh.coord(vertexOld,:);
                coordI   = obj.mesh.coord(vertexI,:);
                coordNew = obj.mesh.coord(vertexNew,:);
                
                uOld = obj.computeUnitVector(coordI,coordOld);
                uNew = obj.computeUnitVector(coordI,coordNew);
                c = obj.computeAllCellsOfVertex(vertexI);
                cellB  = cellBar(c,:);
                uIuBar = obj.computeUnitVector(coordI,cellB);
            
            
                alpha1 = obj.computeAngle(uOld);
                alpha2 = obj.computeAngle(uNew);
                angle  = obj.computeAngle(uIuBar);
                
                angleRot = alpha1;
            
            
            
                angle  = mod(angle  - angleRot,2*pi);
                alpha1 = mod(alpha1 - angleRot,2*pi);
                alpha2 = mod(alpha2 - angleRot,2*pi);
                
                angleN = [];
                isRight = false;
                for i=1:length(angle)
                    aV(i,:) = [cos(angle(i)),sin(angle(i))];
                    angleN(i) = obj.computeAngle(aV(i,:));
                    angleRef = angleN(i);
                    isRight(i,1) = angleRef < alpha2 && angleRef > alpha1;
                end
                
                cellRight = [cellRight; c(isRight)];
                cellLeft  = [cellLeft ;  c(~isRight)];
            
            end
            
            cellRight = unique(cellRight);
            cellLeft = unique(cellLeft);
            
            figure
            hold on
            obj.plotMesh();
            obj.plotStraightPath();
            obj.plotClosestVertex();
            obj.plotVerticesPath();
            
            s.coord = obj.mesh.coord;
            s.connec = obj.mesh.connec(cellRight,:);
            m = Mesh(s);
            m.plot()
            
            
            figure
            hold on
            obj.plotMesh();
            obj.plotStraightPath();
            obj.plotClosestVertex();
            obj.plotVerticesPath();
            
            s.coord = obj.mesh.coord;
            s.connec = obj.mesh.connec(cellLeft,:);
            m = Mesh(s);
            m.plot()            
            
        end
        
        function angle = computeAngle(obj,aV)
            angle = atan2(aV(:,2),aV(:,1));
            angle = mod(angle,2*pi);
        end
        
     
        
        function u = computeUnitVector(obj,coordA,coordB)
            u(:,1) = coordB(:,1) - coordA(:,1);
            u(:,2) = coordB(:,2) - coordA(:,2);
            nU = u(:,1).^2 + u(:,2).^2;
            u(:,1) = u(:,1)./nU;
            u(:,2) = u(:,2)./nU;            
        end
        
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
        
        function computePath(obj)            
            i = 1;
            vertex       = obj.closestVertex;
            pVertexes    = vertex;
            pEdges       = zeros(1);
            pCells       = {obj.computeAllCellsOfVertex(vertex)};
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
                pCells{i}    = obj.computeAllCellsOfVertex(newVertex);                
            end
            obj.pathVertices = pVertexes;
            obj.pathEdges    = pEdges;
            obj.pathCells    = obj.computeUniqueCells(pCells);
        end
        

        function edges = computeAllEdgesOfVertex(obj,vertex)
            vertexInEdges = obj.mesh.edges.nodesInEdges;            
            isInEdge = vertexInEdges == vertex;            
            allEdges  = 1:size(isInEdge,1);
            allEdgesA(:,1) = allEdges(isInEdge(:,1));
            allEdgesB(:,1) = allEdges(isInEdge(:,2));
            edges = [allEdgesA;allEdgesB];            
        end
        
        function cells = computeAllCellsOfVertex(obj,vertex)
            vertexInCell  = obj.mesh.connec;            
            isInCell      = any(vertexInCell == vertex,2);            
            allCells(:,1) = 1:size(isInCell,1);
            cells         = allCells(isInCell);                        
        end
        
        function oV = computeOtherVertexOfEdge(obj,edge,vertex)
            vertexesInEdges = obj.mesh.edges.nodesInEdges;
            vertexesOfEdge  = vertexesInEdges(edge,:);
            oVB = setdiff(vertexesOfEdge(:,2),vertex);
            oVA = setdiff(vertexesOfEdge(:,1),vertex);
            oV = [oVB;oVA];
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
            plot(x,y,'g-','LineWidth',5)
        end
        
    end
    
    methods (Access = private, Static)
        
        function uCells = computeUniqueCells(cells)
            aCells = [];
            for iC = 1:size(cells,2)
                cell = cells{iC};
                nCell = length(cell);
                aCells(end+1:end+nCell) = cell;
            end
            uCells = unique(aCells);            
        end
        
    end
    
    
end