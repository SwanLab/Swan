classdef TopOptComputer < handle

    properties (Access = public)
        computation
        settings
        variables
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = TopOptComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            obj.createSettings();
            topOptSolver = TopOpt_Problem(obj.settings);
            topOptSolver.computeVariables();
            obj.computation = topOptSolver;
            obj.variables.x = topOptSolver.designVariable.fun.fValues;
        end
    end

    methods (Access = private)
          
        function createSettings(obj)
            obj.createOldSettings();
            obj.translateToNewSettings();
        end

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
end