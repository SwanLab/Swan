classdef ReferenceCellOfPathComputer < handle
    
    properties (Access = private)
        isInCell
    end
    
    properties (Access = private)
        pathVertexes
        mesh
    end
    
    methods (Access = public)
        
        function obj = ReferenceCellOfPathComputer(cParams)
            obj.init(cParams)
            
        end
        
        function cellsR = compute(obj)
            cellsR = obj.computeFirstReferenceCell();
            vertex = obj.pathVertexes;        
            for ivertex = 2:length(vertex)
                vI   = vertex(ivertex);
                vOld = vertex(ivertex-1);
                cells    = obj.computeCellsOfVertex(vI);
                cellsOld = obj.computeCellsOfVertex(vOld);
                commonCells = intersect(cells,cellsOld);
                cellsR(ivertex) = commonCells(1);
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.pathVertexes = cParams.pathVertexes;
            obj.mesh         = cParams.mesh;
        end
        
        function cell = computeFirstReferenceCell(obj)
            firstVertex  = obj.pathVertexes(1);         
            cells =  obj.computeCellsOfVertex(firstVertex);
            cell  = cells(1);
        end
        
        function cells = computeCellsOfVertex(obj,vertex)
            vertexInCell  = obj.mesh.connec;
            itIs = (vertexInCell == vertex); 
            cells = any(itIs,2);
            cells = find(cells);
        end            
        
    end
    
end