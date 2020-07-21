classdef CutMeshComputerProvisional < handle
    
    properties (Access = public)
        connec
        coord
        xCoordsIso
        cellContainingSubcell
    end
    
    properties (Access = private)
        cutEdgesComputer  
        cutPointsInElemComputer      
        
        cutEdgesParams
        cutCoordParams
        interiorSubCellsParams
        cutEdgesComputerParams   
        
        cutCoordComputer
    end
    
    properties (Access = private)        
        backgroundMesh
        cutCells
    end
    
    methods (Access = public)
        
        function  obj = CutMeshComputerProvisional(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeCutEdges();
            obj.computeCutCoordinateComputer();  
            obj.coord = obj.cutCoordComputer.coord;
            obj.computeCutPointsInElemComputer();            
            obj.computeConnec();
            obj.computeMesh();
        end        
        
        function m = computeBoundaryMesh(obj)
            s.coord  = obj.cutCoordComputer.xCutPoints;
            s.connec = obj.cutPointsInElemComputer.edgeCutPointInElem;
            s.kFace  = obj.backgroundMesh.kFace -1;
            m = Mesh(s);
        end        
        
        function m = computeMesh(obj)
            sM.connec = obj.connec;
            sM.coord  = obj.coord;
            sM.kFace  = obj.backgroundMesh.kFace;
            m = Mesh(sM);
        end        
        
        function xCutIso = obtainBoundaryXcutIso(obj)
            xCutIso = obj.cutPointsInElemComputer.xCut;
        end
        
        function c = obtainBoundaryCellContainingSubCell(obj)
            c = obj.cutCells;            
        end
        
    end    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.computeAllParams(cParams);
        end
        
        function computeAllParams(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.cutCells       = cParams.cutCells;
            cutCells       = cParams.cutCells;
            ls             = cParams.levelSet;            
            obj.backgroundMesh.computeEdges();
            e = obj.backgroundMesh.edges;
            obj.cutEdgesParams.nodesInEdges = e.nodesInEdges;
            obj.cutEdgesParams.levelSet     = ls;
            obj.cutCoordParams.coord = obj.backgroundMesh.coord;
            obj.cutCoordParams.nodesInEdges = e.nodesInEdges;
            
            cEparams = obj.cutEdgesComputerParams;
            
            cEparams.allNodesinElemParams.finalNodeNumber = size(obj.backgroundMesh.coord,1);
            cEparams.allNodesinElemParams.backgroundConnec = obj.backgroundMesh.connec;
            cEparams.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem;
            cEparams.edgesInElem = e.edgesInElem;
            cEparams.nEdgeByElem = e.nEdgeByElem;
            
            obj.cutEdgesComputerParams = cEparams;
            
            obj.interiorSubCellsParams.isSubCellInteriorParams.levelSet = ls;
            obj.interiorSubCellsParams.cutElems = obj.cutCells;
        end
        
        function computeCutEdges(obj)
            s = obj.cutEdgesParams;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;    
        end
     
        function computeCutCoordinateComputer(obj)
            s = obj.cutCoordParams;
            s.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            cComputer = CutCoordinatesComputer(s);
            cComputer.compute();        
            obj.cutCoordComputer = cComputer;
        end
        
        function computeCutPointsInElemComputer(obj)
            s = obj.cutEdgesComputerParams;
            s.isEdgeCut = obj.cutEdgesComputer.isEdgeCut;
            s.allNodesInElemCoordParams.xCutEdgePoint = obj.cutEdgesComputer.xCutEdgePoint;
            c = CutPointsInElemComputer(s);
            c.compute();
            obj.cutPointsInElemComputer = c;
        end
         
        function computeConnec(obj)
            c = obj.cutPointsInElemComputer;
            sS.bestSubCellCaseSelector.coord = obj.coord;
            sA.subMeshConnecParams = sS;
            sA.xAllNodesInElem = c.xAllNodesInElem;
            sA.allNodesInElem  = c.allNodesInElem;
            sC.isEdgeCutInElem = c.isEdgeCutInElem;
            s = obj.interiorSubCellsParams;          
            sI = s.isSubCellInteriorParams;
            sI.allNodesInElem = c.allNodesInElem;
            s.allSubCellsConnecParams = sA;
            s.subCellsCasesParams = sC; 
            s.isSubCellInteriorParams = sI;
            subCell = InteriorSubCellsConnecComputer(s);
            obj.connec                = subCell.connec;
            obj.xCoordsIso            = subCell.xCoordsIso;
            obj.cellContainingSubcell = subCell.cellContainingSubcell;
        end            

    end
    
end