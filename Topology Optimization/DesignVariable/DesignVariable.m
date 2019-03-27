classdef DesignVariable < handle & matlab.mixin.Copyable
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        value
        meshGiD
    end
    
    methods (Access = public, Abstract)
        
        update(obj,value)
        
    end
    
    methods (Access = public)
        
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
    end
    
end

