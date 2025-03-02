classdef StiffnessTensor2VoigtConverterPS < FourthOrderTensor2VoigtConverterPS
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StiffnessTensor2VoigtConverterPS(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
       
        function selectVoigtTensorClass(obj)
            obj.voigtTensor = StiffnessPlaneStressVoigtTensor();
        end
        
    end
    
    methods (Access = protected,Static)
        
        function f = getVoigtFactor(iv,jv)
           f = 1;
        end

    end
    
    
end
    