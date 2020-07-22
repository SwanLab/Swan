classdef CutMesh < handle
    
    properties (Access = protected)
        backgroundMesh        
        levelSet
        cutCells
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = CutMeshFactory();
            obj = f.create(cParams);            
        end
        
    end
    
    methods (Abstract, Access = public)
        compute(obj)
        
    end
    
    methods (Abstract, Access = protected)
        
        obtainMesh(obj)
        obtainXcoordIso(obj)
        obtainCellContainingSubCells(obj)
        
        obtainBoundaryMesh(obj)                
        obtainBoundaryXcutIso(obj)
        obtainBoundaryCellContainingSubCell(obj)
    end
        
    
    methods (Access = public)
        
        function c = computeInteriorMesh(obj)
            xCoord                   = obj.obtainXcoordIso();
            cellContSubCell          = obj.obtainCellContainingSubCells();            
            c.mesh                   = obj.obtainMesh();
            c.cellContainingSubcell  = cellContSubCell;
            c.cutMeshOfSubCellLocal  = obj.computeCutMeshOfSubCellLocal(xCoord,'INTERIOR');               
        end
        
        function c = computeBoundaryMesh(obj)
            xCoord                   = obj.obtainBoundaryXcutIso(); 
            cellContSubCell          = obj.obtainBoundaryCellContainingSubCell();
            c.mesh                   = obj.obtainBoundaryMesh();
            c.cellContainingSubcell  = cellContSubCell;
            c.cutMeshOfSubCellLocal  = obj.computeCutMeshOfSubCellLocal(xCoord,'BOUNDARY');               
        end
        
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
            %obj.backgroundCutMesh = cParams.backgroundCutMesh;            
            obj.cutCells       = cParams.cutCells;
        end
        
    end
    
  
    
    methods (Access = private, Static)
        
        function m = computeCutMeshOfSubCellLocal(xCoordIso,type)
            coord = xCoordIso;
            nElem = size(coord,3);
            nNode = size(coord,2);
            nDim  = size(coord,1);
            s.coord = reshape(coord,nDim,[])';
            s.connec = reshape(1:nElem*nNode,nNode,nElem)';
            if isequal(type,'INTERIOR')
                kFace = 0;
            else
                kFace = -1;
            end
            s.kFace = kFace;
            m = Mesh(s);
        end
        
    end
    
    
end