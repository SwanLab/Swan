classdef testComputingTopOptWithVademecumData < testShowingError
    
    
    properties (Access = protected)
        tol = 1e-6;
        testName = 'testComputingTopOptWithVademecumData';
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = testComputingTopOptWithVademecumData()
            obj.createTopOptProblem();
        end
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function createTopOptProblem(obj)
            settings = Settings('CantileverTriangleFineM1M2');  
            settingsTopOpt = SettingsTopOptProblem('CaseBenchmark_JSON_B.json',settings);
            
            topOptProblem = TopOpt_Problem(settingsTopOpt);
            topOptProblem.computeVariables();
            topOptProblem.postProcess();
        end
        
        
    end
    
end