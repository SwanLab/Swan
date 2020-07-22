classdef CutMeshProvisionalLine < CutMesh
    
    properties (Access = public)
        connec
        coord
        xCoordsIso
        cellContainingSubcell        
    end
    
    properties (Access = private)
        cutEdgesComputer
        cutCoordComputer
    end
    
    properties (Access = private)
        isoCoord
    end
    
    methods (Access = public)
        
        function obj = CutMeshProvisionalLine(cParams)
            obj.init(cParams)
            obj.isoCoord = [-1 1];
        end
        
        function compute(obj)
            obj.computeCutCoordinate();
            obj.coord = obj.cutCoordComputer.coord;
            obj.computeConnec();
            obj.computeXcoordIso();
            obj.cellContainingSubcell = obj.cutCells;
        end
        
    end
    
    methods (Access = protected)
        
        function m = obtainMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            s.kFace  = obj.backgroundCutMesh.kFace;
            m = Mesh(s);
        end
        
        function x = obtainXcoordIso(obj)
            x = obj.xCoordsIso;
        end        
        
        function c = obtainCellContainingSubCells(obj)
           c = obj.cellContainingSubcell; 
        end
        
        function obtainBoundaryMesh(obj)
          
        end         
        
        function obtainBoundaryXcutIso(obj)
                   
        end
        
        function  obtainBoundaryCellContainingSubCell(obj)
          
        end                  
        
    end
    
    methods (Access = private)

        function computeCutCoordinate(obj)
            s.levelSet      = obj.levelSet;
            s.nodesInEdges  = obj.backgroundCutMesh.connec;  
            obj.cutEdgesComputer = CutEdgesComputer(s);
            obj.cutEdgesComputer.compute();    
            
            s.coord            = obj.backgroundCutMesh.coord;
            s.nodesInEdges     = obj.backgroundCutMesh.connec;
            s.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            cComputer = CutCoordinatesComputer(s);
            cComputer.compute();        
            obj.cutCoordComputer = cComputer;            
        end
        
        function computeConnec(obj)
            nodeA = obj.backgroundCutMesh.connec(:,1);
            nodeB = obj.backgroundCutMesh.connec(:,2);
            lsA = obj.levelSet(nodeA);
            isNodeAinterior = lsA < 0;
            isNodeBinterior = ~isNodeAinterior;
            interiorNode(isNodeAinterior,1) = nodeA(isNodeAinterior);
            interiorNode(isNodeBinterior,1) = nodeB(isNodeBinterior);
            lastNode = size(obj.backgroundCutMesh.coord,1);
            nElem    = obj.backgroundCutMesh.nelem;
            boundaryNode(:,1) = lastNode + (1:nElem);
            obj.connec = [interiorNode boundaryNode];
        end
        
        function computeXcoordIso(obj)
            nodeA = obj.backgroundCutMesh.connec(:,1);
            lsA = obj.levelSet(nodeA);
            isNodeAinterior = lsA < 0;         
            isNodeBinterior = ~isNodeAinterior;            
            nDim = 1;
            nnode = 2;
            nelem = length(obj.cutCells);
            xIso = zeros(nDim,nnode,nelem);
            xIso(1,1,:) = obj.cutEdgesComputer.xCutEdgePoint;
            xIso(1,2,isNodeAinterior) = obj.isoCoord(1);
            xIso(1,2,isNodeBinterior) = obj.isoCoord(2);
            obj.xCoordsIso = xIso;
        end
        
    end
    
end