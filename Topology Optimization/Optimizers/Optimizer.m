classdef Optimizer < handle
    
    properties (Access = protected)
        stop_vars
        hasConverged
    end
    
    properties (Access = public)
        target_parameters = struct;
    end
        
end