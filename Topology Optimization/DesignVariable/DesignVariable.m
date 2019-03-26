classdef DesignVariable < handle
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        value
    end
    
    methods (Access = public, Abstract)
        
        update(obj,value)
        
    end
    
end

