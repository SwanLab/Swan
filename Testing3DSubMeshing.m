classdef Testing3DSubMeshing < handle
    
    properties (Access = public)
    end
    
    properties (Access = private)
        backgroundMesh
        boundaryMesh
        uMesh
        levelSet
    end
    
    properties (Access = private)
        coord        
    end
    
    methods (Access = public)
        
        function obj = Testing3DSubMeshing()
            obj.init();
            obj.createLevelSet();
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
            obj.createUnfittedMesh();
            obj.plotUnfittedMesh();
            obj.computeCutPoints();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.coord = [0 0 0;
                     1 0 0;
                     0 1 0;
                     0 0 1];    
            %obj.coord = rand(4,3);
        end
        
        function createLevelSet(obj)
             b = -10;
             a = 10;
            % obj.levelSet = b + (a-b)*rand(size(obj.coord,1),1);           
             obj.levelSet = [-7.8496;-9.7731;-8.3404;8.3622];
        end
        
        function createBackgroundMesh(obj)
            s.connec = [1 2 3 4];
            s.coord  = obj.coord;
            m = Mesh(s);
            obj.backgroundMesh = m;
        end
        
        function createBoundaryMesh(obj)
            connec = [1 2 3;
                      1 2 4;
                      1 3 4;
                      2 3 4];
            for iFace = 1:size(connec,1)
               con = connec(iFace,:);
               s.coord =  obj.coord(con,:);
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
        end
        
        function computeCutPoints(obj)
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
            
            
            sS.bestSubCellCaseSelector.coord = obj.coord;
            sA.subMeshConnecParams = sS;
            sA.xAllNodesInElem = cE.xAllNodesInElem;
            sA.allNodesInElem  = cE.allNodesInElem;
            
            nodes = obj.backgroundMesh.connec;
            ls = zeros(size(nodes));
            for iNode = 1:size(nodes,2)
                ls(:,iNode) = obj.levelSet(nodes(:,iNode));                 
            end            
            cutCase = 1 - heaviside(ls);            
            
            
            sC.cutCase = cutCase;

            
            cutCells = 1;
            
            
            s.connec = obj.backgroundMesh.connec;
            s.levelSet = obj.levelSet;
            subCells = SubCellsCasesComputer(s);
            subCells.compute();
            subCellCases = subCells.subCellCases;            
            
            
            sI.isSubCellInteriorParams.levelSet = obj.levelSet;
            
             
            sI.subCellCases   = subCellCases;
            sI.allNodesInElem = cE.allNodesInElem;
            s.allSubCellsConnecParams = sA;
            s.subCellsCasesParams = sC; 
            s.isSubCellInteriorParams = sI;
            s.cutElems = cutCells;            
            subCell = InteriorSubCellsConnecComputer(s);
            
            
              
        
            
            
            
        end
        
        
    end
    
end