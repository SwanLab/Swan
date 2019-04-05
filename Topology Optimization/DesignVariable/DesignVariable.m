classdef DesignVariable < handle & matlab.mixin.Copyable
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        value
        meshGiD
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        type
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

