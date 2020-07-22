classdef CutMesh < handle
    
    properties (GetAccess = protected, SetAccess = private)
        backgroundMesh        
        levelSet
        cutCells
    end
    
    properties (Access = protected)
        mesh
        xCoordsIso
        cellContainingSubcell  
        boundaryMesh
        cellContainingSubCellBoundary
        xCoordsIsoBoundary        
    end
    
    properties (SetAccess = protected, GetAccess = public) 
        innerCutMesh
        boundaryCutMesh
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
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.cutCells       = cParams.cutCells;
        end
        
        function computeInnerCutMesh(obj)
            s.mesh                  = obj.mesh;
            s.xCoordsIso            = obj.xCoordsIso;
            s.cellContainingSubcell = obj.cellContainingSubcell;
            obj.innerCutMesh = InnerCutMesh(s);
        end
        
        function computeBoundaryCutMesh(obj)
            s.mesh                  = obj.boundaryMesh;
            s.xCoordsIso            = obj.xCoordsIsoBoundary;
            s.cellContainingSubcell = obj.cellContainingSubCellBoundary;
            obj.boundaryCutMesh = BoundaryCutMesh(s);
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