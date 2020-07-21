classdef CutMeshProvisionalOthers < handle
    
    properties (Access = public)
       coord
       connec
       xCoordsIso
       cellContainingSubcell
       
       type
    end
    
    properties (Access = private)
        subcellsMesher
        cutPointsCalculator        
        memoryManager
        ndimUnf
        
        nCutCells
    end
    
    properties (Access = private)
        backgroundCutCells
        backgroundMesh
        backgroundGeomInterpolation
        levelSet
    end
    
    methods (Access = public)
        
        function obj = CutMeshProvisionalOthers(cParams)
            obj.init(cParams)    
            obj.createSubCellsMesher();
            obj.createMemoryManager();
            obj.createCutPointsCalculator();            
        end
        
        function compute(obj)
            obj.computeCutMesh();
        end
        
        function m = computeMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            if isequal(obj.type,'INTERIOR')
                s.kFace = obj.backgroundMesh.kFace;
            else
                s.kFace = obj.backgroundMesh.kFace -1;
            end
            m = Mesh(s);
        end
        
        function m = computeBoundaryMesh(obj)
            m = obj.computeMesh();
        end
        
        function x = obtainBoundaryXcutIso(obj)
            x = obj.xCoordsIso;
        end
        
        function c = obtainBoundaryCellContainingSubCell(obj)
            c = obj.cellContainingSubcell;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type                        = cParams.type;            
            obj.levelSet                    = cParams.levelSet;            
            obj.backgroundCutCells          = cParams.cutCells;
            obj.backgroundMesh              = cParams.backgroundMesh;
            obj.backgroundGeomInterpolation = cParams.interpolationBackground;   
            obj.nCutCells                   = length(obj.backgroundCutCells);           
        end
        
        function createSubCellsMesher(obj)
            sS.ndimIso            = obj.backgroundMesh.geometryType;
            sS.type               = obj.type;
            sS.posNodes           = obj.backgroundGeomInterpolation.pos_nodes;
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

        function computeCutMesh(obj)
            obj.computeSubcells();
            obj.computeGlobalUnfittedMesh();
            obj.computeXcoordIso();            
            obj.cellContainingSubcell = obj.memoryManager.cellContainingSubcell;            
        end
        
        function computeXcoordIso(obj)
            subCell = obj.memoryManager.subcellIsoCoords;
            obj.xCoordsIso = permute(subCell,[3 2 1]);                        
        end
        
        function computeGlobalUnfittedMesh(obj)
            obj.computeGlobalCoordinates();
            obj.computeGlobalConnectivities();
        end
        
        function obj = computeSubcells(obj)
            obj.memoryManager.allocateMemory();
            obj.computeCutPoints();            
            for icut = 1:obj.nCutCells %Vectorize
                icell = obj.backgroundCutCells(icut);               
                newSubcells = obj.computeThisCellSubcells(icut,icell);                
                newCellContainingNodes   = repmat(icell,[newSubcells.nNodes 1]);
                newCellContainingSubcell = repmat(icell,[newSubcells.nSubcells 1]);                
                obj.memoryManager.saveNewSubcells(newSubcells,newCellContainingNodes,newCellContainingSubcell);
            end
            obj.memoryManager.freeSpareMemory();
            
        end
        
        function computeCutPoints(obj)
            s.meshBackground              = obj.backgroundMesh;
            s.levelSet_background         = obj.levelSet;
            s.backgroundCutCells          = obj.backgroundCutCells;
            s.backgroundGeomInterpolation = obj.backgroundGeomInterpolation;
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
        
        function computeGlobalConnectivities(obj)
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
                indexes = obj.findIndexesComparingCoords(coordsSubCell,obj.coord);
                conn(isub,:) = indexes(connecLocal(isub,:));
            end
            obj.connec = conn;
        end
        
        function computeGlobalCoordinates(obj)
            [coord,ind1,ind2] = unique(obj.memoryManager.coord_global_raw,'rows','stable');
            obj.coord = coord;
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