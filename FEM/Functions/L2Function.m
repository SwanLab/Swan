classdef L2Function < handle
    
    properties (Constant, Access = public)
        fType = 'L2';
    end
    
    properties (Access = protected)
        mesh
    end
    
    methods (Access = public)
        function fun = project(obj,target)
            s.mesh          = obj.mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end
    end

end