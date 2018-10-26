classdef StrainTensor2VoigtConverter < SecondOrderTensor2VoigtConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StrainTensor2VoigtConverter(Tensor)
           obj.computeConversion(Tensor) 
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
    end
    
end

