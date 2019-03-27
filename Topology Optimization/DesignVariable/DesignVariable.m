classdef DesignVariable < handle
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        value
        meshGiD
    end
    
    methods (Access = public, Abstract)
        
        update(obj,value)
        
    end
    
end

