classdef testExplorSettNumHomogDefault < testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    properties (Access = protected)
        tol = 1e-13;
        testName = 'testExplorSettNumHomogDefault';
        variablesToStore = {'Ch'};
    end
    
    properties (Access = private)
        homogenizer
        numHomogParams        
        fileMicroTestName
    end
    
    methods (Access = public)
        
        function obj = testExplorSettNumHomogDefault()
            obj.init();
            obj.computeChomog();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.numHomogParams = SettingsNumericalHomogenizer();
        end
        
        function computeChomog(obj)
            d = obj.numHomogParams;
            obj.homogenizer = NumericalHomogenizer(d);
            obj.computedVar{1} = obj.homogenizer.Ch;            
        end
    end
    
end