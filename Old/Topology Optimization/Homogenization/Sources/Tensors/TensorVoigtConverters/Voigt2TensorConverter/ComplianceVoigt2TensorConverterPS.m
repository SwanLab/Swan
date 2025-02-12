classdef ComplianceVoigt2TensorConverterPS < FourthOrderVoigt2TensorConverterPS
    
    properties
    end
    
    methods (Access = public)
        
        function obj = ComplianceVoigt2TensorConverterPS(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected,Static)
        
        function f = getVoigtFactor(iv,jv)
            if ((iv > 2 && jv <= 2) )  || ((iv <= 2 && jv > 2) )
                f = 2;
            elseif (iv>2 && jv>2)
                f = 4;
            else
                f = 1;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function selectTensorClass(obj)
            obj.tensor = CompliancePlaneStressTensor();
        end
        
    end
    
    
end
    