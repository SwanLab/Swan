classdef CutMesh < Mesh
    
    properties (GetAccess = public, SetAccess = private)
         subcellIsoCoords
         cellContainingSubcell
    end
    
    properties (Access = private)
        backgroundMesh
    end
        
    methods (Access = public)
        
        function obj = CutMesh(cParams)
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.subcellIsoCoords = cParams.subcellIsoCoords;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
            obj.computeDescriptorParams();
        end
        
        function b = getBackgroundMesh(obj)
            b = obj.backgroundMesh;
        end
        
    end      
    
    methods (Access = protected)
        
        function computeEmbeddingDim(obj)
            obj.embeddedDim = obj.ndim; %max(1,obj.ndim - 1);            
        end
        
    end    
    
end