classdef NewProblemSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
    end
    
    properties (Access = private)
        
    end
    
    methods (Static, Access = public)

        function obj = create(s)
            s.problemType = 'QuadraticOptimization';
            s.solverType = 'Direct';
            s.useSchur = true;
            obj = QuadraticOptimizationSolver(s);
        end

    end

    methods (Access = public)
        
        function obj = NewProblemSolver(cParams)
            obj.init(cParams);
        end

        function [u,L] = solve(obj)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type               = cParams.solverType;
            obj.mode               = cParams.solverMode;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.BCApplier          = cParams.BCApplier;
            obj.solver             = cParams.solver;
        end
        
    end
    
end