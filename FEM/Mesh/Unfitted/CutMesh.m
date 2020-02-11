classdef CutMesh < Mesh
    
    properties (GetAccess = public, SetAccess = private)
        subcellIsoCoords
        cellContainingSubcell
        globalConnec
    end
    
    properties (Access = private)
        backgroundMesh
        
        backgroundFullCells
        backgroundEmptyCells
    end
    
    properties (GetAccess = ?CutPointsCalculator_Abstract, SetAccess = private)
        backgroundGeomInterpolation
        backgroundCutCells
        levelSet_background
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
    end
    
    methods (Access = public)
        
        function obj = CutMesh(cParams)
            obj.type   = cParams.unfittedType;
            obj.backgroundMesh = cParams.meshBackground;
            
%             obj.build(cParams);
%             obj.init(cParams);
%             
%             obj.subcellsMesher.link(obj);
%             obj.memoryManager.link(obj);
%             if obj.isLevelSetCrossingZero()
%                 obj.computeCutMesh();
%             else
%                 obj.returnNullMesh();
%             end
%             
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
            obj.subcellIsoCoords = cParams.subcellIsoCoords;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
            
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
            obj.computeGlobalConnec();
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
        
        function build(obj,cParams)
            obj.ndim = cParams.meshBackground.ndim;
            builder = UnfittedMesh_Builder_Factory.create(cParams.unfittedType,obj.ndim);
            builder.build(obj);
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
            for ielem = 1:nelem
                icell  = obj.cellContainingSubcell(ielem);
                nodes  = obj.backgroundMesh.connec(icell,:);
                obj.globalConnec(ielem,:) = nodes;
            end
            
        end
    end
    
end