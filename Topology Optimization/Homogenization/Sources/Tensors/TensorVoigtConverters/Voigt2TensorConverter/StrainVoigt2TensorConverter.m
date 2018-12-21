classdef StrainVoigt2TensorConverter < SecondOrderVoigt2TensorConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StrainVoigt2TensorConverter(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
        
        function factor = computeVoigtFactor(obj)
            iv = obj.voigtIndex;
            if iv >= 1 && iv <= 3
                factor = 1;
            elseif iv <= 6 && iv > 3
                factor = 2;
            else
                error('VoigtFactor should be between 1 and 6')
            end            
        end
        
        function selectTensorClass(obj)
            obj.tensor = Strain3DTensor();
        end
    end
    
end
