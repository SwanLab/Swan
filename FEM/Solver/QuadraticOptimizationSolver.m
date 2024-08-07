classdef QuadraticOptimizationSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = QuadraticOptimizationSolver(cParams)
            obj.init(cParams)
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fields = cParams.fields;
            obj.LHS    = cParams.LHS;
            obj.RHS    = cParams.RHS;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.useSchur   = cParams.useSchur;
            obj.solverType = cParams.solverType;
        end
        
    end
    
end