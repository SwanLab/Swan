classdef ComputingCutMeshVectorized < handle
    
    properties (Access = public)

    end
    
    properties (Access = private)
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
            obj.createUnfittedMesh();
            obj.plotUnfittedMesh();
            error = obj.computeCutPoints();    
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
        
        function error = computeCutPoints(obj)
            obj.backgroundMesh.computeEdges();
            e = obj.backgroundMesh.edges;
            s.nodesInEdges = e.nodesInEdges;
            s.levelSet     = obj.levelSet;       
            c = CutEdgesComputer(s);
            c.compute();
            cutEdgesComputer = c;   
         
            sC.coord = obj.backgroundMesh.coord;
            sC.nodesInEdges = e.nodesInEdges;  
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
            
            
            cEparams.allNodesinElemParams.finalNodeNumber = size(obj.backgroundMesh.coord,1);
            cEparams.allNodesinElemParams.backgroundConnec = obj.backgroundMesh.connec;
            cEparams.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem;
            cEparams.edgesInElem = e.edgesInElem;
            cEparams.nEdgeByElem = e.nEdgeByElem;            
            cEparams.isEdgeCut = cutEdgesComputer.isEdgeCut;
            cEparams.allNodesInElemCoordParams.xCutEdgePoint = cutEdgesComputer.xCutEdgePoint;
            
            cE = CutPointsInElemComputer(cEparams);
            cE.compute();            
            
            

            
            nodes = obj.backgroundMesh.connec;
            ls = zeros(size(nodes));
            for iNode = 1:size(nodes,2)
                ls(:,iNode) = obj.levelSet(nodes(:,iNode));                 
            end            
            cutCase = 1 - heaviside(ls);            
            
            
            sC.cutCase = cutCase;

            
            cutCells(:,1) = 1:size(nodes,1);
            
            
            s.connec = obj.backgroundMesh.connec;
            s.levelSet = obj.levelSet;
            subCells = SubCellsCasesComputer(s);
            subCells.compute();
            subCellCases = subCells.caseInfo{1}.subCellCases;            

            
            sS.bestSubCellCaseSelector.coord = obj.backgroundMesh.coord;            
            sA.subMeshConnecParams = sS;
            sA.xAllNodesInElem = cE.xAllNodesInElem;
            sA.allNodesInElem  = cE.allNodesInElem;
            sA.subCellCases    = subCellCases ;           

            s.allSubCellsConnecParams = sA;
            s.subCellsCasesParams = sC; 
            s.isSubCellInterior = subCells.caseInfo{1}.isSubCellsInterior;
            s.cutElems = cutCells;    
            
            
            
            switch mode(size(cE.allNodesInElem,2))
                case 7
                    s.nSubCellsByElem = 4;
                case 8
                    s.nSubCellsByElem = 6;
            end            
            
            
            subCell = InteriorSubCellsConnecComputer(s);
            
            sM.connec = subCell.connec;
            sM.coord  = cutCoordComputer.coord;
            
            m = Mesh(sM);
            s.mesh                  = m;
            s.xCoordsIso            = subCell.xCoordsIso;
            s.cellContainingSubcell = subCell.cellContainingSubcell;
            innerCutMesh = InnerCutMesh(s);          
            
            quad = Quadrature.set(innerCutMesh.mesh.type);
            quad.computeQuadrature('CONSTANT');

            connecU = obj.uMesh.innerCutMesh.mesh.connec
            connecI = innerCutMesh.mesh.connec

            
            vR = obj.uMesh.innerCutMesh.mesh.computeDvolume(quad);            
            vA = innerCutMesh.mesh.computeDvolume(quad);                       
            
            vAT = zeros(size(vR));
            vAT(1:length(vA)) = vA;
            volums = [vR; vAT]'
            
            error = abs(sum(vA) - sum(vR))
            
        end
        
    end
    
end