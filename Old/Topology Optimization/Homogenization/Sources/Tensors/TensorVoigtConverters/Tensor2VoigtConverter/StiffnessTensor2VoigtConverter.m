classdef StiffnessTensor2VoigtConverter < FourthOrderTensor2VoigtConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StiffnessTensor2VoigtConverter(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
        
        function selectVoigtTensorClass(obj)
            obj.voigtTensor = Stiffness3DVoigtTensor();
        end  
        
    end
    
    methods (Access = protected,Static)
        
        function f = getVoigtFactor(iv,jv)
           f = 1;
        end
    end
    
    
end
    