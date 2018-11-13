classdef testTopOptComputation < handle
    
    properties (Abstract, Access = protected)
        testName
    end
    
    properties (Access = protected)
       topOpt
    end
    
    methods (Access = protected)
        
        function obj = testTopOptComputation()
           obj.computeVariableThroughTopOptSolver()
           obj.selectComputedVar();
        end
    end
    
    methods (Access = protected)
        
        function computeVariableThroughTopOptSolver(obj)
            file2load = strcat('./Input/',obj.testName);
            settings = Settings(file2load);
            topOptSolver = TopOpt_Problem(settings);
            topOptSolver.preProcess;
            topOptSolver.computeVariables;
            obj.topOpt = topOptSolver;
        end
        
    end
    
    methods (Abstract, Access = protected)
        selectComputedVar(obj)
    end
    
end

