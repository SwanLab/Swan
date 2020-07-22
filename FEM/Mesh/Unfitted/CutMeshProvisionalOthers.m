classdef CutMeshProvisionalOthers < CutMesh
    
    properties (Access = public)
       type
    end
    
    properties (Access = private)        
        subcellsMesher
        cutPointsCalculator        
        memoryManager        
        nCutCells
        coord
        connec
    end
    
    methods (Access = public)
        
        function obj = CutMeshProvisionalOthers(cParams)
            obj.init(cParams);
            obj.nCutCells = length(obj.cutCells);                                              
        end
        
        function compute(obj)
            obj.computeThisInnerCutMesh();    
            if ~isequal(obj.backgroundMesh.geometryType,'Line')                 
                obj.computeThisBoundaryCutMesh();
            end
        end
        
        function computeThisBoundaryCutMesh(obj)
            obj.type  = 'BOUNDARY';
            obj.createSubCellsMesher();
            obj.createMemoryManager();
            obj.createCutPointsCalculator();             
            obj.computeSubcells();
            obj.computeCoord();
            obj.computeConnec();            
            obj.computeMesh();    
            obj.xCoordsIso            = obj.computeXcoordIso();            
            obj.cellContainingSubcell = obj.memoryManager.cellContainingSubcell;             
            obj.computeBoundaryMesh();
            obj.computeBoundaryXCoordsIso();
            obj.computeBoundaryCellContainingSubCell();            
            obj.computeBoundaryCutMesh();                        
        end
        
        function computeThisInnerCutMesh(obj)
            obj.type  = 'INTERIOR';
            obj.createSubCellsMesher();
            obj.createMemoryManager();
            obj.createCutPointsCalculator();             
            obj.computeSubcells();
            obj.computeCoord();
            obj.computeConnec();            
            obj.computeMesh();
            obj.xCoordsIso            = obj.computeXcoordIso();            
            obj.cellContainingSubcell = obj.memoryManager.cellContainingSubcell;  
            obj.computeInnerCutMesh();            
        end
        
    end
  
    methods (Access = private)
        
        function createSubCellsMesher(obj)
            inter = Interpolation.create(obj.backgroundMesh,'LINEAR');
            sS.ndimIso            = obj.backgroundMesh.geometryType;
            sS.type               = obj.type;
            sS.posNodes           = inter.pos_nodes;
            sS.levelSetBackground = obj.levelSet;
            sS.coordsBackground   = obj.backgroundMesh.coord;
            obj.subcellsMesher    = SubcellsMesher.create(sS);
        end           
        
         function createMemoryManager(obj)
            s.ndimIso       = obj.backgroundMesh.geometryType;
            s.unfittedType  = obj.type;
            s.nCutCells     = obj.nCutCells;
            s.ndim          = obj.backgroundMesh.ndim;
            obj.memoryManager  = MemoryManager_MeshUnfitted(s);
        end        
        
        function createCutPointsCalculator(obj)
            obj.cutPointsCalculator  = CutPointsCalculator();            
        end       
        
        function x = computeXcoordIso(obj)
            subCell = obj.memoryManager.subcellIsoCoords;
            x = permute(subCell,[3 2 1]);                        
        end
        
        function computeMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            if isequal(obj.type,'INTERIOR')
                s.kFace = obj.backgroundMesh.kFace;
            else
                s.kFace = obj.backgroundMesh.kFace -1;
            end
            obj.mesh = Mesh(s);            
        end
        
        function obj = computeSubcells(obj)
            obj.memoryManager.allocateMemory();
            obj.computeCutPoints();            
            for icut = 1:obj.nCutCells %Vectorize
                icell = obj.cutCells(icut);               
                newSubcells = obj.computeThisCellSubcells(icut,icell);                
                newCellContainingNodes   = repmat(icell,[newSubcells.nNodes 1]);
                newCellContainingSubcell = repmat(icell,[newSubcells.nSubcells 1]);                
                obj.memoryManager.saveNewSubcells(newSubcells,newCellContainingNodes,newCellContainingSubcell);
            end
            obj.memoryManager.freeSpareMemory();
            
        end
        
        function computeCutPoints(obj)
            s.backgroundMesh              = obj.backgroundMesh;
            s.levelSet_background         = obj.levelSet;
            s.backgroundCutCells          = obj.cutCells;
            obj.cutPointsCalculator.init(s);
            obj.cutPointsCalculator.computeCutPoints();
        end
        
        function subcells = computeThisCellSubcells(obj,icut,icell)
            cutPoints = obj.cutPointsCalculator.getThisCellCutPoints(icut);
            conn      = obj.backgroundMesh.connec(icell,:);
            sS.cellConnec = conn;
            sS.cutPoints  = cutPoints;            
            obj.subcellsMesher.computeSubcells(sS);            
            subcells = obj.subcellsMesher.subcells;
        end
        
        function computeConnec(obj)
            nSubcells = size(obj.memoryManager.connec_local,1);
            cellOfSubCell = obj.memoryManager.cellContainingSubcell;
            coordGlobal   = obj.memoryManager.coord_global_raw;
            connecLocal   = obj.memoryManager.connec_local;
            cellContNodes = obj.memoryManager.cellContainingNodes;
            nnode = size(connecLocal,2);
            conn = zeros(nSubcells,nnode);
            for isub = 1:nSubcells %Vectorize !!!
                cell = cellContNodes == cellOfSubCell(isub);
                coordsSubCell = coordGlobal(cell,:);
                indexes = obj.findIndexesComparingCoords(coordsSubCell);
                conn(isub,:) = indexes(connecLocal(isub,:));
            end
            obj.connec = conn;
        end
        
        function computeCoord(obj)
            allCoord = obj.memoryManager.coord_global_raw;
            obj.coord = unique(allCoord,'rows','stable');            
        end
        
        function I = findIndexesComparingCoords(obj,A)
            B = obj.coord;
            I = zeros(1,size(A,1));
            for inode = 1:size(A,1)
                match = true(size(B,1),1);
                for idime = 1:size(A,2)
                    match = match & B(:,idime) == A(inode,idime);
                end
                I(inode) = find(match,1);
            end
        end
        
        function computeBoundaryMesh(obj)
           obj.boundaryMesh = obj.mesh; 
        end
        
        function computeBoundaryXCoordsIso(obj)
            obj.xCoordsIsoBoundary = obj.xCoordsIso;
        end
        
        function computeBoundaryCellContainingSubCell(obj)
            obj.cellContainingSubCellBoundary = obj.cellContainingSubcell;
        end
        
    end
    
    
end