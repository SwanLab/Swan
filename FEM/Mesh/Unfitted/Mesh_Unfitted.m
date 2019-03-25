classdef Mesh_Unfitted < Mesh ...
        & Mesh_Unfitted_Abstract...
        & Mesh_Unfitted_Properties
    
    methods (Access = public)
        
        function obj = Mesh_Unfitted(cParams)
            if nargin == 0
                cParams = SettingsMeshUnfitted();
            end
            obj.build(cParams);
            obj.init(cParams);
        end
        
        function computeMesh(obj,levelSet)
            obj.updateLevelSet(levelSet);
            obj.classifyCells();
            if obj.isLevelSetCuttingMesh()
                obj.computeUnfittedMesh();
            else
                obj.returnNullMesh();
            end
            obj.computeGeometryType();
        end
        
        function m = computeMass(obj)
            integrator = Integrator.create(obj);
            nnodesBackground = size(obj.levelSet_background);
            M2 = integrator.integrateUnfittedMesh(ones(nnodesBackground),obj);
            m = sum(M2);
        end
        
        function plot(obj)
            h = figure;
            obj.add2plot(axes(h));
            light
            axis equal off
            hold off
        end
        
        function add2plot(obj,ax,removedDim,removedDimCoord)
            meshUnfittedCopy = obj.clone();
            if obj.existPatchingInputs(nargin)
                meshUnfittedCopy = obj.meshPlotter.patchRemovedDimension(meshUnfittedCopy,removedDim,removedDimCoord);
            end
            obj.meshPlotter.plot(meshUnfittedCopy,ax);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            mB = cParams.meshBackground;
            iB = cParams.interpolationBackground;
            obj.ndim = mB.ndim;
            obj.meshBackground = mB;
            obj.backgroundGeomInterpolation = iB;
        end
        
        function build(obj,cParams)
            ndim = cParams.meshBackground.ndim;
            builder = UnfittedMesh_Builder_Factory.create(cParams.unfittedType,ndim);
            builder.build(obj);
        end
        
        function computeUnfittedMesh(obj)
            obj.computeSubcells();
            obj.computeGlobalUnfittedMesh();
        end
        
        function computeGlobalUnfittedMesh(obj)
            obj.computeGlobalCoordinates();
            obj.computeGlobalConnectivities();
        end
        
        function returnNullMesh(obj)
            obj.coord = zeros(0,obj.ndim);
            obj.connec = [];
        end
        
        function obj = computeSubcells(obj)
            obj.memoryManager.allocateMemory();
            
            obj.computeCutPoints();
            
            for icut = 1:obj.nCutCells
                icell = obj.backgroundCutCells(icut);
                
                newSubcells = obj.computeThisCellSubcells(icut,icell);
                
                newCellContainingNodes = repmat(icell,[newSubcells.nNodes 1]);
                newCellContainingSubcell = repmat(icell,[newSubcells.nSubcells 1]);
                
                obj.memoryManager.saveNewSubcells(newSubcells,newCellContainingNodes,newCellContainingSubcell);
            end
            obj.memoryManager.freeSpareMemory();
            obj.memoryManager.transferData();
        end
        
        function computeGlobalCoordinates(obj)
            obj.coord = unique(obj.coord_global_raw,'rows','stable');
        end
        
        function computeGlobalConnectivities(obj)
            nSubcells = size(obj.connec_local,1);
            for isub = 1:nSubcells
                icell = obj.cellContainingSubcell(isub);
                indexes = obj.findSubcellNodesIndexes(icell);
                obj.assembleConnecs(isub,indexes);
            end
        end
        
        function subcells = computeThisCellSubcells(obj,icut,icell)
            cutPoints_thisCell = obj.cutPointsCalculator.getThisCellCutPoints(icut);
            connec_thisCell = obj.meshBackground.connec(icell,:);
            
            obj.subcellsMesher.computeSubcells(connec_thisCell,cutPoints_thisCell);
            
            subcells = obj.subcellsMesher.subcells;
        end
        
        function computeCutPoints(obj)
            obj.cutPointsCalculator.init(obj);
            obj.cutPointsCalculator.computeCutPoints();
        end
        
        function classifyCells(obj)
            [F,E,C] = obj.cellsClassifier.classifyCells(obj.levelSet_background,obj.meshBackground.connec);
            obj.backgroundFullCells = F;
            obj.backgroundEmptyCells = E;
            obj.backgroundCutCells = C;
        end
        
        function updateLevelSet(obj,levelSet_background)
            obj.levelSet_background = levelSet_background;
            obj.subcellsMesher.link(obj);
            obj.memoryManager.link(obj);
        end
        
        function itIs = isLevelSetCuttingMesh(obj)
            itIs = ~isempty(obj.backgroundCutCells);
        end
        
        function indexes = findSubcellNodesIndexes(obj,icell)
            thisSubcellCoords = obj.coord_global_raw(obj.cellContainingNodes == icell,:);
            indexes = obj.findIndexesComparingCoords(thisSubcellCoords,obj.coord);
        end
        
        function assembleConnecs(obj,isub,indexes)
            obj.connec(isub,:) = indexes(obj.connec_local(isub,:));
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

