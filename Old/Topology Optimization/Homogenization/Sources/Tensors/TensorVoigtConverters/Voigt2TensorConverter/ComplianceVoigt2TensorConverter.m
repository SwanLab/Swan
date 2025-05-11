classdef ComplianceVoigt2TensorConverter < FourthOrderVoigt2TensorConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = ComplianceVoigt2TensorConverter(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected,Static)
        
        function f = getVoigtFactor(iv,jv)
            if ((iv > 3 && jv <= 3) )  || ((iv <= 3 && jv > 3) )
                f = 2;
            elseif (iv>3 && jv>3)
                f = 4;
            else
                f = 1;
            end
        end
    end
    
    methods (Access = protected)
        
        function selectTensorClass(obj)
            obj.tensor = Compliance3DTensor();
        end
        
    end
    
end
    