classdef StrainTensor < SecondOrderTensor
   
    properties 
        VoigtFactor
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
        
       
        function makeItPlaneStress(obj)
            obj.transformTensor2Voigt();
            TensorPS = PlaneStressTransformer.transform(obj.tensorVoigt);
            obj.tensorVoigtInPlaneStress = TensorPS;
        end
        
        
        function createWithPlaneStressCompatibility(obj,Ch)
          stress = StressTensor();
          stress.makeItPlaneStress();
          obj.tensorVoigt = (Ch.tensorVoigt)\stress.tensorVoigt;
          
          e = StrainVoigtTensor();
          e.setValue(obj.tensorVoigt);
          obj.tensor = Voigt2TensorConverter.convert(e);
          TensorPS = PlaneStressTransformer.transform(obj.tensorVoigt);
          obj.tensorVoigtInPlaneStress = TensorPS;
        end
        
        
    end   
   
end

