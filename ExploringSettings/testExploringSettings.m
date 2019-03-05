classdef testExploringSettings < testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    properties (Access = protected)
        tol = 1e-13;
        testName = 'testExploringSettings';
        variablesToStore = {'Ch'};
    end
    
    properties (Access = private)
        homogenizer
        numHomogDataBase        
        fileMicroTestName
    end
    
    methods (Access = public)
        
        function obj = testExploringSettings()
            obj.init();
            obj.computeNumericalHomogenizerDataBase();
            obj.computeChomog();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileMicroTestName = 'test2d_micro';
        end
        
        function computeNumericalHomogenizerDataBase(obj)
            defaultDB = NumericalHomogenizerDataBase([obj.fileMicroTestName,'.m']);
            dB = defaultDB.dataBase;
            dB.outFileName                   = obj.fileMicroTestName;
            dB.print                         = true;
            dB.levelSetDataBase.levelSetType = 'full';
            obj.numHomogDataBase = dB;
        end
        
        function computeChomog(obj)
            d = obj.numHomogDataBase;
            obj.homogenizer = NumericalHomogenizer(d);
            obj.computedVar{1} = obj.homogenizer.Ch;            
        end
    end
    
end