classdef InnerCutMesh < handle
    
    properties (Access = public)
        mesh
        xCoordsIso
        cellContainingSubcell
        cutMeshOfSubCellLocal
        cutMeshOfSubCellGlobal
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = InnerCutMesh(cParams)
            obj.init(cParams)
            obj.computeCutMeshOfSubCellLocal();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh                  = cParams.mesh;
            obj.xCoordsIso             = cParams.xCoordsIso;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;            
        end
        
        function computeCutMeshOfSubCellLocal(obj)
            coord = obj.xCoordsIso;
            nElem = size(coord,3);
            nNode = size(coord,2);
            nDim  = size(coord,1);
            s.coord = reshape(coord,nDim,[])';
            s.connec = reshape(1:nElem*nNode,nNode,nElem)';
            s.kFace = 0;          
            m = Mesh(s);
            obj.cutMeshOfSubCellLocal = m;
        end   
    
        
        
    end
    
end