classdef testTopOptComputation < handle
    
    properties (Abstract, Access = protected)
        testName
    end
    
    properties (Access = protected)
       topOpt
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
            topOptSolver.computeVariables();
            obj.topOpt = topOptSolver;
        end
          
        function createSettings(obj)
            obj.createOldSettings();
            obj.translateToNewSettings();
        end
        
    end
    
    methods (Access = private)
        
        function createOldSettings(obj)
            fileName = obj.testName;
            s = Settings(fileName);             
            s.warningHoleBC = false;
            s.printIncrementalIter = false; 
            s.printChangingFilter = false;
            s.printing = false;
            obj.settings = s;
        end
        
        function translateToNewSettings(obj)
            translator = SettingsTranslator();
            translator.translate(obj.settings);
            fileName = translator.fileName;
            obj.settings  = SettingsTopOptProblem(fileName);
        end
        
    end
    
    methods (Abstract, Access = protected)
        selectComputedVar(obj)
    end
    
end

