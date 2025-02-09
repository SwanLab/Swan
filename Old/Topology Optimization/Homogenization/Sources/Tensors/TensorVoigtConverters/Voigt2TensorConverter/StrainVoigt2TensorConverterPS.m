classdef StrainVoigt2TensorConverterPS < SecondOrderVoigt2TensorConverterPS
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StrainVoigt2TensorConverterPS(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
        
        function factor = computeVoigtFactor(obj)
            iv = obj.voigtIndex;
            if iv >= 1 && iv <= 2
                factor = 1;
            elseif iv <= 3
                factor = 2;
            else
                error('VoigtFactor should be between 1 and 3')
            end            
        end
        
        function selectTensorClass(obj)
            obj.tensor = StrainPlaneStressTensor();
        end
    end
    
end
