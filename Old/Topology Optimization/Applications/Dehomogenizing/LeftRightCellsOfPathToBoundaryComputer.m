classdef LeftRightCellsOfPathToBoundaryComputer < handle
    
    properties (Access = private)
        allBaricentersCoord
        cellsOfVertex
        isRight
        isCellLeft
        isCellRight
    end

    properties (Access = private)
        mesh
        pathVertexes  
        singularElement
        isCoherent
    end    
    
    methods (Access = public)
        
        function obj = LeftRightCellsOfPathToBoundaryComputer(cParams)
            obj.init(cParams)            
        end
        
        function [cR,cL] = compute(obj)
            if length(obj.pathVertexes) == 1
            else
            obj.computeAllBaricenterCoords();
            obj.computeSingularityCell();           
            obj.computeInitCells();
            obj.computeIntermidiateCells();
            obj.computeFinalCells()
            end
            cR = obj.isCellRight;
            cL = obj.isCellLeft;            
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
            obj.singularElement = cParams.singularElement;
            obj.isCoherent      = cParams.isCoherent;
            obj.isCellRight = false(obj.mesh.nelem,1);            
            obj.isCellLeft  = false(obj.mesh.nelem,1);            
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
            isNotSinVer = obj.cellsOfVertex~=obj.singularElement;
            obj.cellsOfVertex = obj.cellsOfVertex(isNotSinVer);

            aPast = obj.computeSingularityToFirstVertexAngle(v1);            
            aNext = obj.computeEdgeAngle(v1,v2);
            aBar = obj.computeVertexToBaricenterAngle(v1);
            itIs = mod(aBar-aPast,2*pi) < mod(aNext-aPast,2*pi); % & mod(aBar-aPast,2*pi)~=0;            
          %  isNotSingElem =  obj.cellsOfVertex ~= obj.singularElement;            
            obj.isRight = itIs;% & isNotSingElem;
            obj.appendCellsOnRight();
            obj.appendCellsOnLeft();  
       end   

       function computeSingularityCell(obj)
            v1 = obj.pathVertexes(1);    
            v2 = obj.pathVertexes(2);    
            v0 = obj.computeVertex0();

            a10 = obj.computeEdgeAngle(v1,v0);
            a12 = obj.computeEdgeAngle(v1,v2);
            obj.cellsOfVertex = obj.singularElement;
            aBar = obj.computeVertexToBaricenterAngle(v1);            
            itIs =  mod(aBar-a10,2*pi) < mod(a12-a10,2*pi); 


            % vO = obj.computeOtherCoherentVertex();
            % vEnd = obj.pathVertexes(end);
            % obj.cellsOfVertex = obj.singularElement;
            % aO1 = obj.computeEdgeAngle(v1,vO);
            % a1End = obj.computeEdgeAngle(v1,vEnd);
            % aBar = obj.computeVertexToBaricenterAngle(v1);            
            % itIs =  mod(aBar-a1End,2*pi) < mod(aO1-a1End,2*pi); 
            % 
            
            obj.isRight = itIs;
            obj.appendCellsOnRight();
            obj.appendCellsOnLeft();             
       end

       function oV = computeOtherCoherentVertex(obj)
            v1 = obj.pathVertexes(1);
            vS = obj.mesh.connec(obj.singularElement,:);
            vOthers = vS ~= v1;
            isCV = obj.isCoherent.getFvaluesDisc();
            isCS = isCV(1,:,obj.singularElement);            
            isCS = squeeze(isCS);
            isCO = isCS & vOthers;
            oV = vS(isCO);
            oV = oV(end);
       end

       function v0 = computeVertex0(obj)
           isS  = obj.singularElement;
           isCV = obj.isCoherent.getFvaluesByElem();
           isCS = isCV(1,:,isS);
           isCS = squeeze(isCS);
           vertexS = obj.mesh.connec(isS,:);
           if all(isCS)
               for idof = 1:obj.mesh.nnodeElem
                   iVertex = vertexS(idof);
                   cells = obj.mesh.computeAllCellsOfVertex(iVertex);
                   isCo = squeeze(obj.isCoherent.fValues(1,:,cells));
                   isCellC = all(isCo);
                   notCoh(idof) = sum(~isCellC);
               end
               [~,isV0] = max(notCoh);
           else
               v1 = obj.pathVertexes(1);
               isV1 = vertexS == v1;
               isV0 = ~isV1 & ~isCS;
           end
           v0 = vertexS(isV0);
       end

     

       function a = computeSingularityToFirstVertexAngle(obj,v1)
            coordS  = obj.allBaricentersCoord(obj.singularElement,:);
            coord1 = obj.mesh.coord(v1,:);            
            u      = obj.computeUnitVector(coord1,coordS); 
            a  = obj.computeAngle(u);             
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
        
        function computeCellsOfVertex(obj,vertex)
            c = obj.mesh.computeAllCellsOfVertex(vertex);  
            obj.cellsOfVertex = c;
        end        
        
        function appendCellsOnRight(obj)
            obj.isCellRight(obj.cellsOfVertex) = obj.isRight;            
        end
        
        function appendCellsOnLeft(obj)
            isLeft = ~obj.isRight;
            obj.isCellLeft(obj.cellsOfVertex) = isLeft;            
        end
        
        function plotRightCells(obj)
            figure();
            obj.mesh.plot();
            obj.plotVerticesPath();
            obj.plotLineOfCells(obj.isCellRight);                        
        end
        
        function plotLeftCells(obj)
            figure();
            obj.mesh.plot();
            obj.plotVerticesPath();
            obj.plotLineOfCells(obj.isCellLeft);            
        end        
        
        function plotLineOfCells(obj,cells)
            s.coord  = obj.mesh.coord;
            s.connec = obj.mesh.connec(cells,:);
            m = Mesh.create(s);
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