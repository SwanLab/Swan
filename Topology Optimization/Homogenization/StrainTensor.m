classdef StrainTensor < SecondOrderTensor
   
    properties

    end
    
    methods 
        
        function obj = StrainTensor()
          obj = obj@SecondOrderTensor();
        end
        
        
        
        function computeVoigtFactor(obj,iv)
            if iv >= 1 && iv <= 3
                obj.VoigtFactor = 1;
            elseif iv <= 6 && iv > 3
                obj.VoigtFactor = 2;
            else
                error('VoigtFactor should be between 1 and 6')
            end
            
        end
    end   
   
end

