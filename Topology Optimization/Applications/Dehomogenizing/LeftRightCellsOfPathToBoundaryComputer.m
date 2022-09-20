classdef LeftRightCellsOfPathToBoundaryComputer < handle
    
    properties (Access = private)
        mesh
        pathVertexes        
    end
    
    properties (Access = private)
        allBaricentersCoord
        cellsOfVertex
        isRight
        cellRight
        cellLeft
    end
    
    methods (Access = public)
        
        function obj = LeftRightCellsOfPathToBoundaryComputer(cParams)
            obj.init(cParams)            
        end
        
        function [cR,cL] = compute(obj)
           obj.computeAllBaricenterCoords();
            obj.computeInitCells();
            obj.computeIntermidiateCells();
            obj.computeFinalCells()
            obj.makeCellsUnique();
            cR = obj.cellRight;
            cL = obj.cellLeft;
        end        
        
        function plot(obj)
           obj.plotRightCells();
           obj.plotLeftCells();
        end
        
    end
    
    methods (Access = private)        
        
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.pathVertexes = cParams.pathVertexes;
        end     
        
        function computeAllBaricenterCoords(obj)
           bCoord = obj.mesh.computeBaricenter()'; 
           obj.allBaricentersCoord = bCoord;
        end
        
        function computeCellsOnRight(obj,vI,vNext,vPast)
            aPast = obj.computeEdgeAngle(vI,vPast);            
            aNext = obj.computeEdgeAngle(vI,vNext);            
            aBar  = obj.computeVertexToBaricenterAngle(vI);            
            itIs = mod(aBar-aPast,2*pi) < mod(aNext-aPast,2*pi);
            obj.isRight = itIs;
        end    
        
       function computeInitCells(obj)
            vertex = obj.pathVertexes;           
            v1 = vertex(1);
            v2 = vertex(2);
            obj.computeCellsOfVertex(v1);                        
            obj.computeFinalCellsOnRight(v2,v1);
            obj.appendCellsOnRight();
            obj.appendCellsOnLeft();             
       end        
        
       function computeIntermidiateCells(obj)
           vertex = obj.pathVertexes;
           for i = 2:(length(vertex)-1)
               vOld = vertex(i-1);
               v    = vertex(i);
               vNew = vertex(i+1);
               obj.computeCellsOfVertex(v);
               obj.computeCellsOnRight(v,vNew,vOld);
               obj.appendCellsOnRight();
               obj.appendCellsOnLeft();
           end
       end
        
        function computeFinalCells(obj)
            vertex = obj.pathVertexes;            
            vOld   = vertex(end-1);            
            vFinal = vertex(end);            
            obj.computeCellsOfVertex(vFinal);                        
            obj.computeFinalCellsOnRight(vFinal,vOld);
            obj.appendCellsOnRight();
            obj.appendCellsOnLeft();             
        end
        
        function computeFinalCellsOnRight(obj,v,vOld)
            aPast = obj.computeEdgeAngle(v,vOld);            
            aNext = obj.computeEdgeAngle(vOld,v);
            aBar = obj.computeVertexToBaricenterAngle(v);
            itIs = mod(aBar-aPast,2*pi) < mod(aNext-aPast,2*pi);
            obj.isRight = itIs;  
        end              
        
        function alpha = computeEdgeAngle(obj,vertexI,otherVertex)
            otherCoord = obj.mesh.coord(otherVertex,:);
            coordI     = obj.mesh.coord(vertexI,:);
            u          = obj.computeUnitVector(coordI,otherCoord);  
            alpha      = obj.computeAngle(u);            
        end
        
        function alpha = computeVertexToBaricenterAngle(obj,vertexI)
            bCoord = obj.allBaricentersCoord;
            cells  = obj.cellsOfVertex;
            coordB = bCoord(cells,:);            
            coordI = obj.mesh.coord(vertexI,:);
            u      = obj.computeUnitVector(coordI,coordB);           
            alpha = obj.computeAngle(u);                        
        end        
        
        function makeCellsUnique(obj)
            obj.cellRight = unique(obj.cellRight);
            obj.cellLeft = unique(obj.cellLeft);            
        end
        
        function computeCellsOfVertex(obj,vertex)
            c = obj.mesh.computeAllCellsOfVertex(vertex);  
            obj.cellsOfVertex = c;
        end        
        
        function appendCellsOnRight(obj)
            c = obj.cellsOfVertex(obj.isRight);
            n = length(c);            
            obj.cellRight(end+1:end+n) = c;
        end
        
        function appendCellsOnLeft(obj)
            isLeft = ~obj.isRight;
            c = obj.cellsOfVertex(isLeft);
            n = length(c);
            obj.cellLeft(end+1:end+n) = c;
        end
        
        function plotRightCells(obj)
            figure();
            obj.mesh.plot();
            obj.plotVerticesPath();
            obj.plotLineOfCells(obj.cellRight);                        
        end
        
        function plotLeftCells(obj)
            figure();
            obj.mesh.plot();
            obj.plotVerticesPath();
            obj.plotLineOfCells(obj.cellLeft);            
        end        
        
        function plotLineOfCells(obj,cells)
            s.coord  = obj.mesh.coord;
            s.connec = obj.mesh.connec(cells,:);
            m = Mesh(s);
            m.plot()                           
        end
               
        function plotVerticesPath(obj)
            cV = obj.pathVertexes;
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            plot(x,y,'g-','LineWidth',5)
        end        
        
    end
    
    methods (Access = private, Static)
        
        function angle = computeAngle(aV)
            angle = atan2(aV(:,2),aV(:,1));
            angle = mod(angle,2*pi);
        end                
        
        function u = computeUnitVector(coordA,coordB)
            u(:,1) = coordB(:,1) - coordA(:,1);
            u(:,2) = coordB(:,2) - coordA(:,2);
            nU = u(:,1).^2 + u(:,2).^2;
            u(:,1) = u(:,1)./nU;
            u(:,2) = u(:,2)./nU;
        end        
    
        
    end
    
end