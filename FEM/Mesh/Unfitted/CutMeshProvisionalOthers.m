classdef CutMeshProvisionalOthers < CutMesh
    
    properties (Access = public)
       xCoordsIso
       cellContainingSubcell       
       mesh
       
       type
    end
    
    properties (Access = private)        
        subcellsMesher
        cutPointsCalculator        
        memoryManager        
        nCutCells
    end
    
    methods (Access = public)
        
        function obj = CutMeshProvisionalOthers(cParams)
            obj.init(cParams);
            obj.type                        = cParams.type;                        
            obj.nCutCells                   = length(obj.cutCells);                       
            obj.createSubCellsMesher();
            obj.createMemoryManager();
            obj.createCutPointsCalculator();            
        end
        
        function compute(obj)
            obj.computeSubcells();
            obj.mesh                 = obj.computeGlobalUnfittedMesh();
            obj.xCoordsIso           = obj.computeXcoordIso();            
            obj.cellContainingSubcell = obj.memoryManager.cellContainingSubcell;           
        end
        
    end
    
    methods (Access = protected)
        
        function m = obtainMesh(obj)
            m = obj.mesh;
        end
        
        function x = obtainXcoordIso(obj)
            x = obj.xCoordsIso;
        end     
        
        function c = obtainCellContainingSubCells(obj)
           c = obj.cellContainingSubcell; 
        end        
        
        function m = obtainBoundaryMesh(obj)
            m = obj.mesh;
        end
        
        function x = obtainBoundaryXcutIso(obj)
            x = obj.xCoordsIso;
        end
        
        function c = obtainBoundaryCellContainingSubCell(obj)
            c = obj.cellContainingSubcell;
        end
        
    end
    
    methods (Access = private)
        
        function createSubCellsMesher(obj)
            inter = Interpolation.create(obj.backgroundMesh,'LINEAR');
            sS.ndimIso            = obj.backgroundMesh.geometryType;
            sS.type               = obj.type;
            sS.posNodes           = inter.pos_nodes;
            sS.levelSetBackground = obj.levelSet;
            sS.coordsBackground   = obj.backgroundCutMesh.coord;
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
        
        function m = computeGlobalUnfittedMesh(obj)
            s.coord  = obj.computeGlobalCoordinates();
            s.connec = obj.computeGlobalConnectivities(s.coord);
            if isequal(obj.type,'INTERIOR')
                s.kFace = obj.backgroundCutMesh.kFace;
            else
                s.kFace = obj.backgroundCutMesh.kFace -1;
            end
            m = Mesh(s);            
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
        
        function conn = computeGlobalConnectivities(obj,coord)
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
                indexes = obj.findIndexesComparingCoords(coordsSubCell,coord);
                conn(isub,:) = indexes(connecLocal(isub,:));
            end
        end
        
        function coord = computeGlobalCoordinates(obj)
            [coord,ind1,ind2] = unique(obj.memoryManager.coord_global_raw,'rows','stable');            
        end
        
    end
    
    methods (Access = private, Static)
        
        function I = findIndexesComparingCoords(A,B)
            I = zeros(1,size(A,1));
            for inode = 1:size(A,1)
                match = true(size(B,1),1);
                for idime = 1:size(A,2)
                    match = match & B(:,idime) == A(inode,idime);
                end
                I(inode) = find(match,1);
            end
        end
        
    end
    
    
end