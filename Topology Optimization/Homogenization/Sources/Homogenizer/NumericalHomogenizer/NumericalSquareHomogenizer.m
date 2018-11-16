classdef NumericalSquareHomogenizer < NumericalHomogenizer
    
    properties
    end
    
    methods (Access = public)
        
        function obj = NumericalSquareHomogenizer(fileName,print)
            obj.init(fileName,print);
            obj.generateMicroProblem();
            obj.computeHomogenizedVariables();
            obj.createDensityPrinter()
            obj.print()
        end
        
    end
    
    
    methods (Access = protected)
        
        function createDensity(obj)
            
            
            
            
            
            
            obj.densityCreator = DensityCreatorByInitializer(obj.microProblem);
        end
        
    end
    
end

