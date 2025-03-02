classdef STLExporter < handle

    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        
        function obj = STLExporter(cParams)
            obj.init(cParams)
        end

        function export(obj)
            s2g = SwanGiDInterface();
            s2g.exportSTL(obj.mesh);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
        end
        
    end
    
end