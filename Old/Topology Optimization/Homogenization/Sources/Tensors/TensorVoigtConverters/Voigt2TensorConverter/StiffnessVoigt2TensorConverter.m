classdef StiffnessVoigt2TensorConverter < FourthOrderVoigt2TensorConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StiffnessVoigt2TensorConverter(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
        
        function selectTensorClass(obj)
            obj.tensor = Stiffness3DTensor();
        end  
        
    end
    
    methods (Access = protected,Static)
        
        function f = getVoigtFactor(iv,jv)
           f = 1;
        end
    end
    
    
end
    