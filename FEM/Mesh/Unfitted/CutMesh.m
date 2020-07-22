classdef CutMesh < handle
    
    properties (GetAccess = protected, SetAccess = private)
        backgroundMesh        
        levelSet
        cutCells
    end
    
    properties (SetAccess = protected, GetAccess = public)
        mesh
        xCoordsIso
        cellContainingSubcell          
    end
    
    properties (GetAccess = protected, SetAccess = protected)        
        boundaryMesh
        cellContainingSubCellBoundary
        xCoordsIsoBoundary              
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
    
    methods (Access = public)
        
        function c = computeInteriorMesh(obj)
            xCoord                   = obj.xCoordsIso;
            cellContSubCell          = obj.cellContainingSubcell;            
            c.mesh                   = obj.mesh();
            c.cellContainingSubcell  = cellContSubCell;
            c.cutMeshOfSubCellLocal  = obj.computeCutMeshOfSubCellLocal(xCoord,'INTERIOR');               
        end
        
        function c = computeBoundaryMesh2(obj)
            xCoord                   = obj.xCoordsIsoBoundary; 
            cellContSubCell          = obj.cellContainingSubCellBoundary;
            c.mesh                   = obj.boundaryMesh;
            c.cellContainingSubcell  = cellContSubCell;
            c.cutMeshOfSubCellLocal  = obj.computeCutMeshOfSubCellLocal(xCoord,'BOUNDARY');               
        end
        
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
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