classdef testExplorSettNumHomogCustom < testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    properties (Access = protected)
        tol = 1e-13;
        testName = 'testExplorSettNumHomogCustom';
        variablesToStore = {'Ch'};
    end
    
    properties (Access = private)
        homogenizer
        numHomogParams        
        fileMicroTestName
    end
    
    methods (Access = public)
        
        function obj = testExplorSettNumHomogCustom()
            obj.init();
            obj.computeChomog();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.numHomogParams = SettingsNumericalHomogenizer('paramsNumericalHomogenizer_example');
        end
        
        function computeChomog(obj)
            d = obj.numHomogParams;
            obj.homogenizer = NumericalHomogenizer(d);
            obj.computedVar{1} = obj.homogenizer.Ch;            
        end
    end
    
end