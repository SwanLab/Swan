classdef StressVoigt2TensorConverter < SecondOrderVoigt2TensorConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StressVoigt2TensorConverter(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
        
        function factor = computeVoigtFactor(obj)
            factor = 1;               
        end            
        
        
        function selectTensorClass(obj)
            obj.tensor = Stress3DTensor();
        end
    end
    
end
