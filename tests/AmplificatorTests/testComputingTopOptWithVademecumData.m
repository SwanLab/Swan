classdef testComputingTopOptWithVademecumData < testShowingError ...
   
    
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
            fileName = 'CantileverTriangleCoarse_Case_1_1_1';
            settings = Settings(fileName);

            settings.vademecumFileName = 'SmoothRectangle';
            settings.homegenizedVariablesComputer = 'ByVademecum';%'ByVademecum';
            settings.materialInterpolation = 'SIMP-ALL';
            settings.designVariable = 'MicroParams';'Density';'MicroParams';
            
            settings.optimizer = 'PROJECTED GRADIENT';

            topOptProblem = TopOpt_Problem(settings);
            
            
        end
   
        
    end
    
end