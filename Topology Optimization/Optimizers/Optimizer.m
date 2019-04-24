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
    
    methods (Access = public, Static)
       
        function obj = create(cParams)
            f = OptimizerFactory();
            obj = f.create(cParams);            
        end
        
    end
        
end