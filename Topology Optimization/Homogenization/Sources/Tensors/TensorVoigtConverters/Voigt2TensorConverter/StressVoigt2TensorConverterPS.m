classdef StressVoigt2TensorConverterPS < SecondOrderVoigt2TensorConverterPS
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StressVoigt2TensorConverterPS(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
        
        function factor = computeVoigtFactor(obj)
            factor = 1;               
        end            
        
        function selectTensorClass(obj)
            obj.tensor = StressPlaneStressTensor();
        end
    end
    
end
