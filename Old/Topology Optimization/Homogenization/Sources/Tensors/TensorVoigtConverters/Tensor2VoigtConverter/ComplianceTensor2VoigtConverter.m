classdef ComplianceTensor2VoigtConverter < FourthOrderTensor2VoigtConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = ComplianceTensor2VoigtConverter(tensor)
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
        
        function selectVoigtTensorClass(obj)
            obj.voigtTensor = Compliance3DVoigtTensor();
        end
        
    end
    
    
end