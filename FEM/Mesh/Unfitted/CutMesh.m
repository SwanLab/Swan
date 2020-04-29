classdef CutMesh < Mesh
    
    properties (GetAccess = public, SetAccess = private)
        cellContainingSubcell
    end
    
    properties (Access = private)
        subcellIsoCoords
        
        
        backgroundMesh
        
        backgroundEmptyCells
        backgroundCutCells
        backgroundFullCells
        
        
        levelSet_unfitted
        
        index1
        index2
        
        subcellsMesher
        cutPointsCalculator
        meshPlotter
        cellsClassifier
        memoryManager
        
        coord_iso
        coord_global_raw
        connec_local
        cellContainingNodes
        
        ndimUnf
        
        nCutCells
        ndimBackground
        
        levelSet_background
        backgroundGeomInterpolation
    end
    
    methods (Access = public)
        
        function obj = CutMesh(cParams)
            obj.meshBackground = cParams.meshBackground;
            obj.isInBoundary   = cParams.isInBoundary;
            obj.type           = cParams.type;
            obj.ndim           = cParams.meshBackground.ndim;
            
            isTriangle = isequal(cParams.meshBackground.geometryType,'TRIANGLE');
            isInterior = isequal(obj.type,'INTERIOR');            
            cutElems = cParams.backgroundCutCells;                   
            thereIsCutElem = ~isempty(cutElems);
            isQuad = isequal(cParams.meshBackground.geometryType,'QUAD');
            if isTriangle && isInterior && thereIsCutElem 
                 
                ls = cParams.levelSet;
                connecCut = cParams.meshBackground.connec(cutElems,:);            
    
                coord = cParams.meshBackground.coord(:,1:2);     
                
                s.coord  = coord;
                s.connec = connecCut;
                backgroundCutMesh = Mesh().create(s);
                %backgroundCutMesh.computeEdges();
                
                s.backgroundMesh = backgroundCutMesh;
                s.cutElems = cutElems;
                s.levelSet = ls;
                
                cM = CutMeshComputerProvisional(s);
                
                
                
                cM.compute();
               
                %obj.coord = zeros(size(cParams.meshBackground.coord));
                obj.connec = cM.connec;
                obj.coord = zeros(size(cM.coord,1),size(cParams.meshBackground.coord,2));
                obj.coord(:,1:2)  = cM.coord;
                obj.subcellIsoCoords = permute(cM.xCoordsIso,[1 3 2]);
                obj.cellContainingSubcell = cM.cellContainingSubcell;
                
            elseif 0%isQuad && isInterior && thereIsCutElem 
                
                ls = cParams.levelSet;
                connecCut = cParams.meshBackground.connec(cutElems,:);            
    
                coord = cParams.meshBackground.coord(:,1:2);     
                
                s.coord  = coord;
                s.connec = connecCut;
                backgroundCutMesh = Mesh().create(s);
                %backgroundCutMesh.computeEdges();
                
                s.backgroundMesh = backgroundCutMesh;
                s.cutElems = cutElems;
                s.levelSet = ls;
                s.lastNode = max(cParams.meshBackground.connec(:));

                cM = CutMeshProvisionalQuadrilater(s);                
                cM.compute();
               
                %obj.coord = zeros(size(cParams.meshBackground.coord));
                obj.connec = cM.connec;
                obj.coord = zeros(size(cM.coord,1),size(cParams.meshBackground.coord,2));
                obj.coord(:,1:2)  = cM.coord;
                obj.subcellIsoCoords = permute(cM.xCoordsIso,[1 3 2]);
                obj.cellContainingSubcell = cM.cellContainingSubcell;                
                
                
                
            else
                
        
                obj.init(cParams);
                obj.build();
                
                
                
                if obj.isLevelSetCrossingZero()
                    obj.computeCutMesh();
                else
                    obj.returnNullMesh();
                end
                
                
            end
            
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
            
        end
        
        
        function setLevelSetUnfitted(obj,LS)
            obj.levelSet_unfitted = LS;
        end
        
        function add2plot(obj,ax,removedDim,removedDimCoord)
            meshUnfittedCopy = obj.clone();
            if obj.existPatchingInputs(nargin)
                meshUnfittedCopy = obj.meshPlotter.patchRemovedDimension(meshUnfittedCopy,removedDim,removedDimCoord);
            end
            bF = obj.backgroundFullCells;
            obj.meshPlotter.plot(meshUnfittedCopy,ax,bF);
        end
        
        function xGauss = computeIsoGaussPoints(obj,quad)
%             coord = obj.subcellIsoCoords;
%             coord = permute(coord,[1 3 2]);                        
%             shape = obj.createShapes(quad.posgp);
%             nDime = size(coord,2);
%             nNode = obj.nnode;
%             nElem = size(coord,1);
%             nGaus = quad.ngaus;
%             xGauss = zeros(nGaus,nElem,nDime);
%             for kNode = 1:nNode
%                 shapeKJ(:,1) = shape(kNode,:);
%                 xKJ(1,:,:) = coord(:,:,kNode);
%                 xG = bsxfun(@times,shapeKJ,xKJ);
%                 xGauss = xGauss + xG;
%             end
%             xGauss = permute(xGauss,[3 1 2]);
%             
            coord = obj.subcellIsoCoords;
            coord = permute(coord,[3 2 1]);                                           
            int = Interpolation.create(obj,'LINEAR');
            int.computeShapeDeriv(quad.posgp);
            shape = int.shape;            
            nDime = size(coord,1);
            nNode = obj.nnode;
            nElem = size(coord,3);
            nGaus = quad.ngaus;
            xGauss = zeros(nDime,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shape(kNode,:);
                xKJ = coord(:,kNode,:);
                xG = bsxfun(@times,shapeKJ,xKJ);
                xGauss = xGauss + xG;
            end
        end
        
        function shapes = createShapes(obj,xG)
            int = Interpolation.create(obj,'LINEAR');
            int.computeShapeDeriv(xG);
            shapes = permute(int.shape,[1 3 2]);
        end        
        
    end
    
    methods (Access = protected)
        
        function computeEmbeddingDim(obj)
            switch obj.type
                case 'BOUNDARY'
                    obj.embeddedDim = obj.ndim - 1;
                case {'INTERIOR'}
                    obj.embeddedDim = obj.ndim;
                otherwise
                    error('EmbeddedDim not defined')
            end
            
            if obj.isInBoundary
                obj.embeddedDim = obj.ndim - 1;
            end
        end
        
    end
    
    methods (Access = private)
        
        function createNdimUnf(obj)
            if obj.isInBoundary
                nUnf = obj.ndim - 1;
            else
                nUnf = obj.ndim;
            end
            obj.ndimUnf = nUnf;
        end
        
        function init(obj,cParams)
            obj.levelSet_background   = cParams.levelSet;
            obj.backgroundFullCells   = cParams.backgroundFullCells;
            obj.backgroundEmptyCells  = cParams.backgroundEmptyCells;
            obj.backgroundCutCells    = cParams.backgroundCutCells;
            obj.nCutCells             = length(obj.backgroundCutCells);
            
            obj.backgroundMesh              = cParams.meshBackground;
            obj.backgroundGeomInterpolation = cParams.interpolationBackground;
        end
        
        function createSubCellsMesher(obj)
            sS.ndimIso = obj.ndimUnf;
            sS.type = obj.type;
            
            sS.posNodes = obj.backgroundGeomInterpolation.pos_nodes;
            sS.levelSetBackground = obj.levelSet_background;
            sS.coordsBackground = obj.meshBackground.coord;
            obj.subcellsMesher = SubcellsMesher.create(sS);
        end
        
        
        
        function createMeshPlotter(obj)
            s.ndimIso  = obj.ndimUnf;
            s.type     = obj.type;
            obj.meshPlotter = MeshPlotter.create(s);
        end
        
        function createMemoryManager(obj)
            s.ndimIso       = obj.ndimUnf;
            s.unfittedType  = obj.type;
            s.nCutCells     = obj.nCutCells;
            s.ndim          = obj.ndim;
            obj.memoryManager  = MemoryManager_MeshUnfitted(s);
        end
        
        function build(obj)
            obj.createNdimUnf();
            obj.createSubCellsMesher();
            obj.createMeshPlotter();
            obj.createMemoryManager();
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
                
                newCellContainingNodes   = repmat(icell,[newSubcells.nNodes 1]);
                newCellContainingSubcell = repmat(icell,[newSubcells.nSubcells 1]);
                
                nnode = (newSubcells.nNodes);
                
                
                obj.memoryManager.saveNewSubcells(newSubcells,newCellContainingNodes,newCellContainingSubcell);
            end
            obj.memoryManager.freeSpareMemory();
            
            
            obj.coord_iso             = obj.memoryManager.coord_iso;
            obj.coord_global_raw      = obj.memoryManager.coord_global_raw;
            obj.subcellIsoCoords      = obj.memoryManager.subcellIsoCoords;
            obj.connec_local          = obj.memoryManager.connec_local;
            obj.connec                = obj.memoryManager.connec;
            obj.levelSet_unfitted     = obj.memoryManager.levelSet_unfitted;
            obj.cellContainingNodes   = obj.memoryManager.cellContainingNodes;
            obj.cellContainingSubcell = obj.memoryManager.cellContainingSubcell;
            
        end
        
        function computeCutPoints(obj)
            s.meshBackground              = obj.meshBackground;
            s.levelSet_background         = obj.levelSet_background;
            s.backgroundCutCells          = obj.backgroundCutCells;
            s.backgroundGeomInterpolation = obj.backgroundGeomInterpolation;
            obj.cutPointsCalculator.init(s);
            obj.cutPointsCalculator.computeCutPoints();
        end
        
        function subcells = computeThisCellSubcells(obj,icut,icell)
            cutPoints_thisCell = obj.cutPointsCalculator.getThisCellCutPoints(icut);
            connec_thisCell = obj.meshBackground.connec(icell,:);
            
            
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