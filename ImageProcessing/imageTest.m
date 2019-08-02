classdef test < handle
    
    properties (Access = private)
        savedData
        computedData
    end
    
    methods (Access = public)
        
        function obj = test(testName)
            obj.loadSavedData(testName);
            obj.computeTest();
            obj.checkEquivalence();
        end
        
    end
    
    methods (Access = private)
        
        function loadSavedData(obj,testName)
            testData = load(testName);
            obj.savedData = testData;
        end
        
        function computeTest(obj,testName)
            problemName  = testName;
            sLoader      = SettingsLoader(problemName);
            settings     = sLoader.settings;            
            imageProblem = DenoisingProblem(settings);
            imageProblem.solve();
            obj.computedData = imageProblem.optimizedImage;
        end
        
    end
    
    
    
    
    
    
end