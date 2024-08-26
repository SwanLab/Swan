classdef QuadraticOptimizationSolver < NewProblemSolver
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        fields
        LHS, RHS
        solverStrategy, solverType
        BCApplier%, boundaryConditions
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = QuadraticOptimizationSolver(cParams)
            obj.init(cParams)
            C = obj.getBoundaryConditionsMatrix();
            canSchur = obj.evaluateSchurUsage();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fields = cParams.fields;
            obj.LHS    = cParams.LHS;
            obj.RHS    = cParams.RHS;
            obj.BCApplier      = cParams.BCApplier;
            obj.solverType     = cParams.solverType;
            obj.solverStrategy = cParams.solverStrategy;
%             obj.boundaryConditions = cParams.boundaryConditions;
        end

        function C = getBoundaryConditionsMatrix(obj)
            C = obj.BCApplier.computeLinearConditionsMatrix('Dirac');
        end

        function canSchur = evaluateSchurUsage(obj)
            canSchur = true;
        end
        
    end
    
end