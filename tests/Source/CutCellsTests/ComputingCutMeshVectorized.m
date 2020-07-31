classdef ComputingCutMeshVectorized < handle
    
    properties (Access = public)

    end
    
    properties (Access = private)
        t2
        t1
        boundaryMesh
        uMesh        
    end
    
    properties (Access = private)
        backgroundMesh
        levelSet
        boundaryConnec
    end
    
    methods (Access = public)
        
        function obj = ComputingCutMeshVectorized(cParams)
            obj.init(cParams)            
        end
        
        function error = compute(obj)
            obj.createBoundaryMesh();
           tic            
            obj.createUnfittedMesh();
          obj.t1 = toc;
           % obj.plotUnfittedMesh();
          tic  
            error = obj.computeCutPoints();    
          obj.t2 =   toc;
          ratio = obj.t1/obj.t2
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.levelSet = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.boundaryConnec = cParams.boundaryConnec;
        end
        
         function createBoundaryMesh(obj)
            connec = obj.boundaryConnec;
            for iFace = 1:size(connec,1)
               con = connec(iFace,:);
               s.coord =  obj.backgroundMesh.coord(con,:);
               s.nodesInBoxFaces = con;
               s.connec = [1 2 3];
               s.kFace = 0;
               s.dimension = [];
               s.isRectangularBox = false;
               m{iFace} = BoundaryMesh(s);
            end
            obj.boundaryMesh = m;
        end
        
        function createUnfittedMesh(obj)
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(obj.levelSet)            
        end
        
        function plotUnfittedMesh(obj)
            obj.uMesh.plotBoundary()
            view([1 1 1])
        end
        
        function subCellCases = computeSubCellCases(obj)
            s.connec     = obj.backgroundMesh.connec;
            s.levelSet   = obj.levelSet;
            subCellCases = SubCellsCasesComputer(s);
            subCellCases.compute();
        end
        
        function error = computeCutPoints(obj)
            
            obj.backgroundMesh.computeEdges();
            
            e = obj.backgroundMesh.edges;
            s.nodesInEdges = e.nodesInEdges;
            s.levelSet     = obj.levelSet;       
            c = CutEdgesComputer(s);
            c.compute();
            cutEdgesComputer = c;   
            
            
            sC.coord            = obj.backgroundMesh.coord;
            sC.nodesInEdges     = e.nodesInEdges;  
            sC.xCutEdgePoint    = cutEdgesComputer.xCutEdgePoint;
            sC.isEdgeCut        = cutEdgesComputer.isEdgeCut;
            cComputer = CutCoordinatesComputer(sC);
            cComputer.compute();        
            cutCoordComputer = cComputer;            
            
            x = cutCoordComputer.xCutPoints(:,1);
            y = cutCoordComputer.xCutPoints(:,2);
            z = cutCoordComputer.xCutPoints(:,3);
            hold on
            plot3(x,y,z,'k*','LineWidth',10,'MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k')
                 
            
            s.edgesInElem   = e.edgesInElem;
            s.isEdgeCut     = cutEdgesComputer.isEdgeCut;
            isEdgeCut       = EdgeCutInElemComputer(s);
            isEdgeCutInElem = isEdgeCut.compute();                   
            
            
            
            subCellCases = obj.computeSubCellCases();
            
            nEdgesCutCase   = [3 4];
            nSubCellsByElem = [4 6];
            
            subCell = cell(length(nEdgesCutCase),1);
            nCutEdges = sum(isEdgeCutInElem,1);
            for icases = 1:length(nEdgesCutCase)
                t = nCutEdges == nEdgesCutCase(icases);
                isEdgeCutInElemCase = isEdgeCutInElem(:,t);
                
                s.isEdgeCutInElem = isEdgeCutInElemCase;
                all2Cut = AllEdges2CutEdgesComputer(s);
                
                cEparams.all2Cut = all2Cut;
                cEparams.allNodesinElemParams.finalNodeNumber = size(obj.backgroundMesh.coord,1);
                cEparams.allNodesinElemParams.connec = obj.backgroundMesh.connec(t,:);
                cEparams.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem(t,:,:);
                cEparams.edgesInElem = e.edgesInElem(t,:);
                cEparams.nEdgeByElem = e.nEdgeByElem;
                cEparams.isEdgeCut = cutEdgesComputer.isEdgeCut;
                cEparams.allNodesInElemCoordParams.xCutEdgePoint = cutEdgesComputer.xCutEdgePoint;
                cEparams.isEdgeCutInElem = isEdgeCutInElemCase;
                cE = CutPointsInElemComputer(cEparams);
                cE.compute();
                cT{icases} = cE;
                
                caseInfo = subCellCases.caseInfo{icases};
                
            
                nodes = obj.backgroundMesh.connec;
                cutCells(:,1) = 1:size(nodes,1);
                
                
                
                sS.bestSubCellCaseSelector.coord = obj.backgroundMesh.coord;
                sA.subMeshConnecParams           = sS;
                sA.xAllNodesInElem               = cE.xAllNodesInElem;
                sA.allNodesInElem                = cE.allNodesInElem;
                sA.subCellCases                  = caseInfo.subCellCases(t,:);
                
                sI.allSubCellsConnecParams = sA;
                sI.isSubCellInterior = caseInfo.isSubCellsInterior(:,t);
                sI.cutElems = cutCells;
                
                sI.nSubCellsByElem = nSubCellsByElem(icases);
                
                if ~isempty(sI.isSubCellInterior)
                subCell{icases} = InteriorSubCellsConnecComputer(sI);
                end
            
            end
            
            connecT = [];
            xCoordsIso = [];
            cellC = [];
            
%
            for icase = 1:2
               subC = subCell{icase};
               if ~isempty(subC)                 
                 connecT = cat(1,connecT,subC.connec);
                 xCoordsIso = cat(3,xCoordsIso,subC.xCoordsIso);

                 cellC = cat(1,cellC,subC.cellContainingSubcell);                 
               end
            end
            
            connec                = connecT;
            coord                 = cutCoordComputer.coord;
            xCoordsIso            = xCoordsIso;
            cellContainingSubcell = cellC;
            
            error = obj.createInnerCutAndPlot(connec,coord,xCoordsIso,cellContainingSubcell);
            
            

            
        end
        
        
        
        function error = createInnerCutAndPlot(obj,connec,coord,xCoordsIso,cellContainingSubcell)            
          
            if ~isempty(obj.uMesh.innerCutMesh)
            quad = Quadrature.set(obj.uMesh.innerCutMesh.mesh.type);
            quad.computeQuadrature('CONSTANT');          
            
            vR = obj.uMesh.innerCutMesh.mesh.computeDvolume(quad);

            sM.connec = connec
            sM.coord  = coord
            
            m = Mesh(sM);
            s.mesh                  = m;
            s.xCoordsIso            = xCoordsIso;
            s.cellContainingSubcell = cellContainingSubcell;
            innerCutMesh = InnerCutMesh(s);

            
            connecU = obj.uMesh.innerCutMesh.mesh.connec
            connecI = innerCutMesh.mesh.connec
            
            
            vA = innerCutMesh.mesh.computeDvolume(quad);
            
            vAT = zeros(size(vR));
            vAT(1:length(vA)) = vA;
            volums = [vR; vAT]'
            
            error = abs(sum(vA) - sum(vR))
            else
                error = 0
            end
        end
        
    end
    
end