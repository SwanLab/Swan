classdef StiffnessVoigt2TensorConverterPS < FourthOrderVoigt2TensorConverterPS
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StiffnessVoigt2TensorConverterPS(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
       
        function selectTensorClass(obj)
            obj.tensor = StiffnessPlaneStressTensor();
        end
        
    end
    
    methods (Access = protected,Static)
        
        function f = getVoigtFactor(iv,jv)
           f = 1;
        end

    end
    
    
end
    