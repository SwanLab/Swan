classdef testComputingFemWithVademecumData < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testComputingFemWithVademecumData';
    end
    
    properties (Access = private)
        fileName
        designVar
        Ctensor
        density
    end
    
    
    methods (Access = public)
        
        function obj = testComputingFemWithVademecumData()
            obj.init();
            obj.createDesignVariableFromRandMxMy();
            obj.computeConstitutiveAndDensityFromVademecum();
        end
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName =  'VademecumSmoothCorner';
        end
        
        function createDesignVariableFromRandMxMy(obj)
            a = 0.01;
            b = 0.99;
            obj.designVar = (b-a).*rand(1000,2) + a;
        end
        
        function computeConstitutiveAndDensityFromVademecum(obj)
            s.fileName = obj.fileName;
            v = VademecumVariablesLoader(s);
            v.load();
            obj.Ctensor = v.Ctensor.compute(obj.designVar);
            obj.density = v.density.compute(obj.designVar);
        end
        
        
    end
    
end