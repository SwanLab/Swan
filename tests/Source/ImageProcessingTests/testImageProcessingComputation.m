classdef testImageProcessingComputation < handle
    
    properties (Abstract, Access = protected)
        testName
    end
    
    properties (Access = protected)
        imageProblem
    end
    
    methods (Access = public)
        
        function obj = testImageProcessingComputation()
            obj.solveImageProblem()
            obj.selectComputedVar();
        end
        
    end
    
    methods (Access = private)
        
        function solveImageProblem(obj)
            problemName  = obj.testName;
            sLoader      = SettingsLoader(problemName);
            settings     = sLoader.settings;
            imProblem = DenoisingProblem(settings);
            imProblem.solve();
            obj.imageProblem = imProblem;
        end
        
    end
    
    methods (Abstract, Access = protected)
        selectComputedVar(obj)
    end
    
end