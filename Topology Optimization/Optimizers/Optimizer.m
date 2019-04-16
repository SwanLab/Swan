classdef Optimizer < handle
    
    properties (Access = protected)
        stop_vars
        has_converged
        constraint_case
        nconstr        
    end
    
    properties (Access = public)
        target_parameters = struct;
    end
    
    methods (Access = public)
        
        function obj = Optimizer(settings)
            obj.nconstr           = settings.nconstr;
            obj.target_parameters = settings.target_parameters;
            obj.constraint_case   = settings.constraint_case;
            obj.has_converged     = false;            
        end
        
    end
    
end