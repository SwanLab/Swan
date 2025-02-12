classdef CutMeshProvisionalLine < CutMesh
    
    properties (Access = private)
        connec
        coord
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
            obj.coord = obj.cutCoordComputer.allValues;
            obj.computeConnec();
            obj.computeXcoordIso();
            obj.cellContainingSubcell = obj.cutCells;
            obj.computeMesh();
            obj.computeInnerCutMesh();
        end
        
    end
    
    methods (Access = private)

        function computeCutCoordinate(obj)
            s.levelSet      = obj.levelSet;
            s.nodesInEdges  = obj.backgroundMesh.connec;
            obj.cutEdgesComputer = CutEdgesComputer(s);
            obj.cutEdgesComputer.compute();
            
            s.fValues       = obj.backgroundMesh.coord;
            s.nodesInEdges     = obj.backgroundMesh.connec;
            s.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            cComputer = CutFunctionValuesComputer(s);
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
            nelem = length(obj.cutCells);
            xIso = zeros(nDim,nnode,nelem);
            xIso(1,1,:) = obj.cutEdgesComputer.xCutEdgePoint;
            xIso(1,2,isNodeAinterior) = obj.isoCoord(1);
            xIso(1,2,isNodeBinterior) = obj.isoCoord(2);
            obj.xCoordsIso = xIso;
        end
        
        function computeMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            s.kFace  = obj.backgroundMesh.kFace;
            obj.mesh = Mesh.create(s);
        end
        
    end
    
end