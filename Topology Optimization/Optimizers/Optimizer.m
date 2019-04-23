classdef Optimizer < handle
    
    properties (Access = public)
        convergenceVars
    end
    
    properties (Access = protected)
        hasConverged
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        name
    end
    
    properties (Access = public)
        target_parameters = struct;
    end
        
end