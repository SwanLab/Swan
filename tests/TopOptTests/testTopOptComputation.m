classdef testTopOptComputation < handle
    
    properties (Abstract, Access = protected)
        testName
    end
    
    properties (Access = protected)
       topOpt
    end
    
    properties (Access = private)
       settings 
    end
    
    methods (Access = protected)
        
        function obj = testTopOptComputation()
           obj.createSettings();
           obj.computeVariableThroughTopOptSolver();
           obj.selectComputedVar();
        end
    end
    
    methods (Access = protected)
        
        function computeVariableThroughTopOptSolver(obj)
            topOptSolver = TopOpt_Problem(obj.settings);
            topOptSolver.preProcess;
            topOptSolver.computeVariables;
            obj.topOpt = topOptSolver;
        end
        
    end
    
    methods (Access = private)
        
        function createSettings(obj)
            file2load = strcat('./Input/',obj.testName);
            sett = Settings(file2load);
            sett.warningHoleBC = false;
            sett.printIncrementalIter = false; 
            sett.printChangingFilter = false;
            obj.settings = sett;
        end
        
    end
    
    methods (Abstract, Access = protected)
        selectComputedVar(obj)
    end
    
end

