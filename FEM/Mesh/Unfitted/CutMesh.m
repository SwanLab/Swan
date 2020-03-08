classdef CutMesh < Mesh
    
    properties (GetAccess = public, SetAccess = private)
        globalConnec
    end
    
    properties (GetAccess = public, SetAccess = ?MemoryManager_MeshUnfitted)
        subcellIsoCoords
        cellContainingSubcell
    end
    
    properties (GetAccess = private, SetAccess = ?MemoryManager_MeshUnfitted)
        coord_iso
        coord_global_raw
        connec_local
        cellContainingNodes
    end
    
    properties (GetAccess = ?PatchedMeshPlotter_Abstract )
        backgroundFullCells
    end
    
    properties (Access = private)
        backgroundMesh
        
        backgroundEmptyCells
        
        levelSet_unfitted
        
        index1
        index2
        
    end
    
    properties (GetAccess = ?CutPointsCalculator, SetAccess = private)
        backgroundCutCells
    end
    
    properties (GetAccess = {?CutPointsCalculator,?SubcellsMesher}, SetAccess = private)
        levelSet_background
        backgroundGeomInterpolation
    end
    
    properties (GetAccess = private, SetAccess = ?UnfittedMesh_AbstractBuilder)
        subcellsMesher
        cutPointsCalculator
        meshPlotter
        cellsClassifier
        memoryManager
    end
    
    properties (GetAccess = ?MemoryManager, SetAccess = ?UnfittedMesh_AbstractBuilder)
        maxSubcells
        nnodesSubcell
    end
    
    properties (GetAccess = ?MemoryManager, SetAccess = private)
        nCutCells
        ndimBackground
        isInBoundary
        
    end
    
    methods (Access = public)
        
        function obj = CutMesh(cParams)
            
            if isfield(cParams,'isInBoundary')
                obj.isInBoundary = cParams.isInBoundary;
            else
                obj.isInBoundary = false;
            end
            
            obj.type   = cParams.unfittedType;
            obj.meshBackground = cParams.meshBackground;
            
            m = obj.meshBackground;
            
            cutConnec = m.connec(cParams.backgroundCutCells,:);
            s.connec = cutConnec;
            nodes = unique(cutConnec(:));
            s.coord = m.coord(nodes,:);
            
            
            
            %             m2 = Mesh().create(s);
            %             figure
            %             m2.plot();
            %             subMesher = SubMesher();
            %             subMesh = subMesher.computeSubMesh(m2);
            %             figure
            %             subMesh.plot();
            
            %ComputeLevelSet
            
            %ComputeMeshUnfitted
            
            %ChangeCparams cParams!!
            
            ndim = cParams.meshBackground.ndim;
            uType = cParams.unfittedType;
            
            builderType = uType;
            if obj.isInBoundary && ndim == 2
               builderType = 'BOUNDARY';
            end
            obj.unfittedType = builderType;            
            
            
            obj.build(ndim);
            obj.init(cParams);
            
            obj.subcellsMesher.link(obj);
            obj.memoryManager.link(obj);
            if obj.isLevelSetCrossingZero()
                obj.computeCutMesh();
            else
                obj.returnNullMesh();
            end
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
            obj.computeGlobalConnec();
        end
        
        
        function setLevelSetUnfitted(obj,LS)
            obj.levelSet_unfitted = LS;
        end
        
        function add2plot(obj,ax,removedDim,removedDimCoord)
            meshUnfittedCopy = obj.clone();
            if obj.existPatchingInputs(nargin)
                meshUnfittedCopy = obj.meshPlotter.patchRemovedDimension(meshUnfittedCopy,removedDim,removedDimCoord);
            end
            obj.meshPlotter.plot(meshUnfittedCopy,ax);
        end
        
    end
    
    methods (Access = protected)
        
        function computeEmbeddingDim(obj)
            switch obj.type
                case 'BOUNDARY'
                    obj.embeddedDim = obj.ndim - 1;
                case {'INTERIOR','COMPOSITE'}
                    obj.embeddedDim = obj.ndim;
                otherwise
                    error('EmbeddedDim not defined')
            end
            isNdim = obj.ndim == 3 || obj.ndim == 2;
            if isequal(obj.type,'INTERIOR') && isNdim && obj.isInBoundary
                obj.embeddedDim = obj.ndim - 1;
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.levelSet_background = cParams.levelSet;
            obj.backgroundFullCells = cParams.backgroundFullCells;
            obj.backgroundEmptyCells = cParams.backgroundEmptyCells;
            obj.backgroundCutCells = cParams.backgroundCutCells;
            obj.nCutCells = length(obj.backgroundCutCells);
            
            obj.backgroundMesh = cParams.meshBackground;
            obj.backgroundGeomInterpolation = cParams.interpolationBackground;
        end
        
        function build(obj,ndim)
            
            obj.ndim = ndim;
            if obj.isInBoundary
                ndimUnf = obj.ndim - 1;
            else
                ndimUnf = obj.ndim;
            end
            
            s.ndimIso = ndimUnf;
            s.unfittedType = obj.unfittedType;
            
            obj.subcellsMesher = SubcellsMesher.create(s);
            obj.meshPlotter    = MeshPlotter.create(s);
            obj.memoryManager  = MemoryManager_MeshUnfitted(s);
            
            obj.cutPointsCalculator  = CutPointsCalculator;
            obj.cellsClassifier      = CellsClassifier;
            
        end
        
        function itIs = isLevelSetCrossingZero(obj)
            itIs = ~isempty(obj.backgroundCutCells);
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
                
                newCellContainingNodes = repmat(icell,[newSubcells.nNodes 1]);
                newCellContainingSubcell = repmat(icell,[newSubcells.nSubcells 1]);
                
                obj.memoryManager.saveNewSubcells(newSubcells,newCellContainingNodes,newCellContainingSubcell);
            end
            obj.memoryManager.freeSpareMemory();
            obj.memoryManager.transferData();
        end
        
        function computeCutPoints(obj)
            obj.cutPointsCalculator.init(obj);
            obj.cutPointsCalculator.computeCutPoints();
        end
        
        function subcells = computeThisCellSubcells(obj,icut,icell)
            cutPoints_thisCell = obj.cutPointsCalculator.getThisCellCutPoints(icut);
            connec_thisCell = obj.meshBackground.connec(icell,:);
            
            obj.subcellsMesher.computeSubcells(connec_thisCell,cutPoints_thisCell);
            
            subcells = obj.subcellsMesher.subcells;
        end
        
        function computeGlobalConnec(obj)
            nnode = obj.backgroundMesh.nnode;
            nelem = obj.nelem;
            obj.globalConnec = zeros(nelem,nnode);
            for ielem = 1:nelem %Vectorize!
                icell  = obj.cellContainingSubcell(ielem);
                nodes  = obj.backgroundMesh.connec(icell,:);
                obj.globalConnec(ielem,:) = nodes;
            end
            
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
                
                %                 p = reshape(obj.index2,[],obj.nCutCells)';
                %                 cell = obj.backgroundCutCells == cellOfSubCell(isub);
                %                 indexes = p(cell,:);
                
                
                %norm(indexes(:) - indexes2(:))
                
                connecV(isub,:) = indexes(connecLocal(isub,:));
                
                
                %  icell = obj.cellContainingSubcell(isub);
                %  indexes = obj.findSubcellNodesIndexes(icell);
                %  obj.assembleConnecs(isub,indexes);
            end
            obj.connec = connecV;
        end
        
        function computeGlobalCoordinates(obj)
            [coord,ind1,ind2] = unique(obj.coord_global_raw,'rows','stable');
            obj.coord = coord;
            obj.index1 = ind1;
            obj.index2 = ind2;
        end
        
        function indexes = findSubcellNodesIndexes(obj,icell)
            thisSubcellCoords = obj.coord_global_raw(obj.cellContainingNodes == icell,:);
            indexes = obj.findIndexesComparingCoords(thisSubcellCoords,obj.coord);
        end
        
        function assembleConnecs(obj,isub,indexes)
            obj.connec(isub,:) = indexes(obj.connec_local(isub,:));
        end
        
        function returnNullMesh(obj)
            obj.coord = zeros(0,obj.ndim);
            obj.connec = [];
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
        
        function theyDo = existPatchingInputs(nInputs)
            if nInputs == 4
                theyDo = true;
            else
                theyDo = false;
            end
        end
        
    end
    
end