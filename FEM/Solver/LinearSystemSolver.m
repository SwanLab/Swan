classdef LinearSystemSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LinearSystemSolver(cParams)
            obj.init(cParams)
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fields = cParams.fields;
            obj.LHS    = cParams.LHS;
            obj.RHS    = cParams.RHS;
            obj.boundaryConditions = cParams.boundaryConditions;
        end
        
    end
    
end