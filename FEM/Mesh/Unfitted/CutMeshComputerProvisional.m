classdef CutMeshComputerProvisional < CutMesh
  
    properties (Access = private)
        connec
        coord        
        
        cutEdgesComputer  
        cutPointsInElemComputer      
        
        cutEdgesParams
        cutCoordParams
        interiorSubCellsParams
        cutEdgesComputerParams   
        
        cutCoordComputer
        
        subCellCases
    end
        
    methods (Access = public)
        
        function  obj = CutMeshComputerProvisional(cParams)
            obj.init(cParams);
            obj.computeAllParams();
        end
        
        function compute(obj)
            obj.computeSubCellCases();            
            obj.computeCutEdges();
            obj.computeCutCoordinateComputer();  
            obj.coord = obj.cutCoordComputer.coord;
            obj.computeCutPointsInElemComputer();            
            obj.computeConnec();
            obj.computeMesh();
            obj.computeBoundaryXCoordsIso();
            obj.computeBoundaryCellContainingSubCell();
            obj.computeBoundaryMesh();                        
            obj.computeInnerCutMesh();
            obj.computeBoundaryCutMesh();            
        end    
        
    end
    
    methods (Access = private)
                
        function computeAllParams(obj)          
            obj.backgroundMesh.computeEdges();
            e = obj.backgroundMesh.edges;
            obj.cutEdgesParams.nodesInEdges = e.nodesInEdges;
            obj.cutEdgesParams.levelSet     = obj.levelSet;
           
            obj.cutCoordParams.coord = obj.backgroundMesh.coord;
            obj.cutCoordParams.nodesInEdges = e.nodesInEdges;
            
            cEparams = obj.cutEdgesComputerParams;
            
            cEparams.allNodesinElemParams.finalNodeNumber = size(obj.backgroundMesh.coord,1);
            cEparams.allNodesinElemParams.backgroundConnec = obj.backgroundMesh.connec;
            cEparams.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem;
            cEparams.edgesInElem = e.edgesInElem;
            cEparams.nEdgeByElem = e.nEdgeByElem;
            
            obj.cutEdgesComputerParams = cEparams;
            
            obj.interiorSubCellsParams.isSubCellInteriorParams.levelSet = obj.levelSet;
            obj.interiorSubCellsParams.cutElems = obj.cutCells;
        end
        
        function computeCutEdges(obj)
            s = obj.cutEdgesParams;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;    
        end
        
        function computeSubCellCases(obj)
            s.connec = obj.backgroundMesh.connec;
            s.levelSet = obj.levelSet;
            subCells = SubCellsCasesComputer(s);
            subCells.compute();
            obj.subCellCases = subCells.subCellCases;
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
            sA.subCellCases    = obj.subCellCases;
            s = obj.interiorSubCellsParams;          
            sI = s.isSubCellInteriorParams;
            sI.allNodesInElem = c.allNodesInElem;
            sI.subCellCases   = obj.subCellCases;
            s.allSubCellsConnecParams = sA;
            s.subCellsCases = obj.subCellCases; 
            s.isSubCellInteriorParams = sI;
            subCell = InteriorSubCellsConnecComputer(s);
            obj.connec                = subCell.connec;
            obj.xCoordsIso            = subCell.xCoordsIso;
            obj.cellContainingSubcell = subCell.cellContainingSubcell;
        end            
        
        function computeMesh(obj)
            sM.connec = obj.connec;
            sM.coord  = obj.coord;
            sM.kFace  = obj.backgroundMesh.kFace;
            obj.mesh = Mesh(sM);            
        end
        
        function computeBoundaryMesh(obj)
            s.coord  = obj.cutCoordComputer.xCutPoints;
            s.connec = obj.cutPointsInElemComputer.edgeCutPointInElem;
            s.kFace  = obj.backgroundMesh.kFace -1;
            obj.boundaryMesh = Mesh(s);
        end                
        
        function computeBoundaryXCoordsIso(obj)
            obj.xCoordsIsoBoundary = obj.cutPointsInElemComputer.xCut;
        end        

        function computeBoundaryCellContainingSubCell(obj)
            obj.cellContainingSubCellBoundary = obj.cutCells;            
        end
    end
    
end