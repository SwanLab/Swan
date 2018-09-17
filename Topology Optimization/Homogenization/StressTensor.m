classdef StressTensor < SecondOrderTensor
   
    properties

    end
    
    methods 
        
        function obj = StressTensor()
          obj = obj@SecondOrderTensor();
        end
        
        
        
        function computeVoigtFactor(obj,iv)
            obj.VoigtFactor = 1;
        end
        
        
        

        
    end
    
end

