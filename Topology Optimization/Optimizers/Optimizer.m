classdef Optimizer < handle
    
    properties (Access = protected)
        stop_vars
        hasConverged
        constraint_case
        nconstr        
    end
    
    properties (Access = public)
        target_parameters = struct;
    end
        
end