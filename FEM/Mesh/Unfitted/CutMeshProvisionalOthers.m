classdef CutMeshProvisionalOthers < handle
    
    properties (Access = public)
       coord
       connec
       xCoordsIso
       cellContainingSubcell
       
       type
    end
    
    properties (Access = private)
        
        coord_iso
        coord_global_raw
        connec_local
        cellContainingNodes        
        
        levelSet_unfitted        
        index1
        index2
        
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
        end
        
        function compute(obj)
            obj.build();            
            obj.computeCutMesh();
            obj.xCoordsIso = permute(obj.xCoordsIso,[3 2 1]);
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
        
        function createNdimUnf(obj)
            obj.ndimUnf = obj.backgroundMesh.geometryType;
        end      
        
         function createMemoryManager(obj)
            s.ndimIso       = obj.ndimUnf;
            s.unfittedType  = obj.type;
            s.nCutCells     = obj.nCutCells;
            s.ndim          = obj.backgroundMesh.ndim;
            obj.memoryManager  = MemoryManager_MeshUnfitted(s);
        end        
        
        function build(obj)
            obj.createNdimUnf();
            obj.createSubCellsMesher();
            obj.createMemoryManager();
            obj.cutPointsCalculator  = CutPointsCalculator;
        end       

        function computeCutMesh(obj)
            obj.computeSubcells();
            obj.computeGlobalUnfittedMesh();
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
                
                nnode = (newSubcells.nNodes);
                obj.memoryManager.saveNewSubcells(newSubcells,newCellContainingNodes,newCellContainingSubcell);
            end
            obj.memoryManager.freeSpareMemory();
            
            
            obj.coord_iso             = obj.memoryManager.coord_iso;
            obj.coord_global_raw      = obj.memoryManager.coord_global_raw;
            obj.xCoordsIso           = obj.memoryManager.subcellIsoCoords;
            obj.connec_local          = obj.memoryManager.connec_local;
            obj.connec                = obj.memoryManager.connec;
            obj.levelSet_unfitted     = obj.memoryManager.levelSet_unfitted;
            obj.cellContainingNodes   = obj.memoryManager.cellContainingNodes;
            obj.cellContainingSubcell = obj.memoryManager.cellContainingSubcell;
            
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
            cutPoints_thisCell = obj.cutPointsCalculator.getThisCellCutPoints(icut);
            connec_thisCell = obj.backgroundMesh.connec(icell,:);
            
            
            sS.cellConnec = connec_thisCell;
            sS.cutPoints = cutPoints_thisCell;
            
            obj.subcellsMesher.computeSubcells(sS);
            
            subcells = obj.subcellsMesher.subcells;
        end
        
        function computeGlobalConnectivities(obj)
            nSubcells = size(obj.connec_local,1);
            cellOfSubCell = obj.cellContainingSubcell;
            coordGlobal   = obj.coord_global_raw;
            connecLocal   = obj.connec_local;
            cellContNodes = obj.cellContainingNodes;
            nnode = size(connecLocal,2);
            connecV = zeros(nSubcells,nnode);
            
            
            for isub = 1:nSubcells %Vectorize !!!
                cell = cellContNodes == cellOfSubCell(isub);
                % size(find(cell))
                coordsSubCell = coordGlobal(cell,:);
                indexes = obj.findIndexesComparingCoords(coordsSubCell,obj.coord);
                
                connecV(isub,:) = indexes(connecLocal(isub,:));
                
                
            end
            obj.connec = connecV;
        end
        
        function computeGlobalCoordinates(obj)
            [coord,ind1,ind2] = unique(obj.coord_global_raw,'rows','stable');
            obj.coord = coord;
            obj.index1 = ind1;
            obj.index2 = ind2;
        end

        
        function createSubCellsMesher(obj)
            sS.ndimIso = obj.ndimUnf;
            sS.type = obj.type;
            sS.posNodes           = obj.backgroundGeomInterpolation.pos_nodes;
            sS.levelSetBackground = obj.levelSet;
            sS.coordsBackground   = obj.backgroundMesh.coord;
            obj.subcellsMesher    = SubcellsMesher.create(sS);
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