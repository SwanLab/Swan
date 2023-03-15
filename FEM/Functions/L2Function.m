classdef L2Function < handle
    
    properties (Constant, Access = public)
        fType = 'L2';
    end
    
    properties (Access = protected)
        mesh
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = L2Function(cParams)
            
        end

        function fun = project(obj,target)
            s.mesh          = obj.mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end
        
    end
    
    methods (Access = private)
        
    end
    
end