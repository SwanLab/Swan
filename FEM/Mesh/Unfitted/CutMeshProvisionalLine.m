classdef CutMeshProvisionalLine < handle
    
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
        backgroundMesh
        cutElems
        levelSet
        isoCoord
    end
    
    methods (Access = public)
        
        function obj = CutMeshProvisionalLine(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeCutCoordinate();
            obj.coord = obj.cutCoordComputer.coord;
            obj.computeConnec();
            obj.computeXcoordIso();
            obj.cellContainingSubcell = obj.cutElems;
        end
        
        function m = computeMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            s.kFace  = obj.backgroundMesh.kFace;
            m = Mesh(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.cutElems       = cParams.cutElems;
            obj.levelSet       = cParams.levelSet;     
            obj.isoCoord       = [-1 1];
        end
        
        function computeCutCoordinate(obj)
            s.levelSet      = obj.levelSet;
            s.nodesInEdges  = obj.backgroundMesh.connec;  
            obj.cutEdgesComputer = CutEdgesComputer(s);
            obj.cutEdgesComputer.compute();    
            
            s.coord = obj.backgroundMesh.coord;
            s.nodesInEdges = obj.backgroundMesh.connec;
            s.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            cComputer = CutCoordinatesComputer(s);
            cComputer.compute();        
            obj.cutCoordComputer = cComputer;            
        end
        
        function computeConnec(obj)
            nodeA = obj.backgroundMesh.connec(:,1);
            nodeB = obj.backgroundMesh.connec(:,2);
            lsA = obj.levelSet(nodeA);
            isNodeAinterior = lsA < 0;
            isNodeBinterior = ~isNodeAinterior;
            interiorNode(isNodeAinterior,1) = nodeA(isNodeAinterior);
            interiorNode(isNodeBinterior,1) = nodeB(isNodeBinterior);
            lastNode = size(obj.backgroundMesh.coord,1);
            nElem    = obj.backgroundMesh.nelem;
            boundaryNode(:,1) = lastNode + (1:nElem);
            obj.connec = [interiorNode boundaryNode];
        end
        
        function computeXcoordIso(obj)
            nodeA = obj.backgroundMesh.connec(:,1);
            lsA = obj.levelSet(nodeA);
            isNodeAinterior = lsA < 0;         
            isNodeBinterior = ~isNodeAinterior;            
            nDim = 1;
            nnode = 2;
            nelem = length(obj.cutElems);
            xIso = zeros(nDim,nnode,nelem);
            xIso(1,1,:) = obj.cutEdgesComputer.xCutEdgePoint;
            xIso(1,2,isNodeAinterior) = obj.isoCoord(1);
            xIso(1,2,isNodeBinterior) = obj.isoCoord(2);
            obj.xCoordsIso = xIso;
        end
        
        
    end
    
end