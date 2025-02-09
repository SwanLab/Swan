classdef InnerCutMesh < handle
    
    properties (Access = public)
        mesh
        xCoordsIso
        cellContainingSubcell
    end
    
   
    methods (Access = public)
        
        function obj = InnerCutMesh(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh                  = cParams.mesh;
            obj.xCoordsIso            = cParams.xCoordsIso;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
        end
            
    end
    
end