classdef ComputingInnerAndBoundaryCutMesh < handle
    
    properties (Access = public)
        innerCutMesh
        cutCoordComputer
    end
    
    properties (Access = private)        
        subCellCases
        cutEdgesComputer
        coord
    end
    
    properties (Access = private)
        backgroundMesh
        levelSet
        boundaryConnec        
    end
    
    methods (Access = public)
        
        function obj = ComputingInnerAndBoundaryCutMesh(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeInnerCutMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.levelSet = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.boundaryConnec = cParams.boundaryConnec;
        end
        
        function computeSubCellCases(obj)
            s.connec     = obj.backgroundMesh.connec;
            s.levelSet   = obj.levelSet;
            obj.subCellCases = SubCellsCasesComputer(s);
            obj.subCellCases.compute();
        end
        
        function computeCutEdges(obj)
            obj.backgroundMesh.computeEdges();
            e = obj.backgroundMesh.edges;
            s.nodesInEdges = e.nodesInEdges;
            s.levelSet     = obj.levelSet;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;
        end
        
        function computeCutCoordinateComputer(obj)
            e = obj.backgroundMesh.edges;
            sC.fValues          = obj.backgroundMesh.coord;
            sC.nodesInEdges     = e.nodesInEdges;
            sC.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            sC.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            cComputer = CutFunctionValuesComputer(sC);
            cComputer.compute();
            obj.cutCoordComputer = cComputer;
        end
        
        function computeInnerCutMesh(obj)
            obj.computeSubCellCases();
            obj.computeCutEdges();
            obj.computeCutCoordinateComputer();
            obj.coord = obj.cutCoordComputer.allValues;
            
            e = obj.backgroundMesh.edges;
            s.edgesInElem   = e.edgesInElem;
            s.isEdgeCut     = obj.cutEdgesComputer.isEdgeCut;
            isEdgeCut       = EdgeCutInElemComputer(s);
            isEdgeCutInElem = isEdgeCut.compute();
            
            e = obj.backgroundMesh.edges;
            nEdgesCutCase   = [2 3 4];
            nSubCellsByElem = [3 4 6];
            
            subCell = cell(length(nEdgesCutCase),1);
            cN = cell(length(nEdgesCutCase),1);
            
            nCutEdges = sum(isEdgeCutInElem,1);
            for icases = 1:length(nEdgesCutCase)
                t = nCutEdges == nEdgesCutCase(icases);
                isEdgeCutInElemCase = isEdgeCutInElem(:,t);
                
                s.isEdgeCutInElem = isEdgeCutInElemCase;
                all2Cut = AllEdges2CutEdgesComputer(s);
                
                cEp.all2Cut = all2Cut;
                cEp.allNodesinElemParams.finalNodeNumber = size(obj.backgroundMesh.coord,1);
                cEp.allNodesinElemParams.connec = obj.backgroundMesh.connec(t,:);
                cEp.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem(t,:,:);
                cEp.edgesInElem = e.edgesInElem(t,:);
                cEp.nEdgeByElem = e.nEdgeByElem;
                cEp.isEdgeCut = obj.cutEdgesComputer.isEdgeCut;
                cEp.allNodesInElemCoordParams.xCutEdgePoint = obj.cutEdgesComputer.xCutEdgePoint;
                cEp.isEdgeCutInElem = isEdgeCutInElemCase;
                cE = CutPointsInElemComputer(cEp);
                cE.compute();
                
                if sum(t) ~= 0
                    
                    cN{icases} = cE;
                    
                    caseInfo = obj.subCellCases.caseInfo{icases};
                    
                    
                    nodes = obj.backgroundMesh.connec;
                    cutCells(:,1) = 1:size(nodes,1);
                    
                    
                    
                    sS.bestSubCellCaseSelector.coord = obj.coord;
                    sA.subMeshConnecParams           = sS;
                    sA.xAllNodesInElem               = cE.xAllNodesInElem;
                    sA.allNodesInElem                = cE.allNodesInElem;
                    sA.subCellCases                  = caseInfo.subCellCases(t,:);
                    
                    sI.allSubCellsConnecParams = sA;
                    sI.isSubCellInterior = caseInfo.isSubCellsInterior(:,t);
                    sI.cutElems = cutCells;
                    
                    sI.nSubCellsByElem = nSubCellsByElem(icases);
                    
                    
                    subCell{icases} = InteriorSubCellsConnecComputer(sI);
                end
                
            end
            
            connecT = [];
            xCoordsIso = [];
            cellC = [];
            
            xCoordsIsoBoundary = [];
            cellContainingSubCellBoundary = [];
            
            for icase = 1:3
                subC = subCell{icase};
                scN  = cN{icases};
                if ~isempty(subC)
                    connecT = cat(1,connecT,subC.connec);
                    xCoordsIso = cat(3,xCoordsIso,subC.xCoordsIso);
                    cellC = cat(1,cellC,subC.cellContainingSubcell);
                end
                
                if ~isempty(scN)
                    xCoordsIsoBoundary = cat(3,xCoordsIsoBoundary,scN.xCutInElem);
                    % cellContainingSubCellBoundary = cat(1,cellContainingSubCellBoundary,obj.cutCells);
                end
                
            end
            
            connec                = connecT;
            % xCoordsIso            = xCoordsIso;
            cellContainingSubcell = cellC;
            
            
            sM.connec = connec;
            sM.coord  = obj.coord;
            
            m = Mesh.create(sM);
            s.mesh                  = m;
            s.xCoordsIso            = xCoordsIso;
            s.cellContainingSubcell = cellContainingSubcell;
            inCutMesh = InnerCutMesh(s);
            inCutMesh.mesh = inCutMesh.mesh.computeCanonicalMesh();
            
            obj.innerCutMesh = inCutMesh;
            
            
        end
        
    end
    
    
    
end