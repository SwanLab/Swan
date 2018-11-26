classdef PST4VoigtFourthOrderTensor < PlaneStressTransformer
    
    
    methods (Access = protected)        
        
        function createPlaneStressTensor(obj)
            obj.psTensor = SymmetricFourthOrderPlaneStressVoigtTensor();
        end
        
    end
   
    
end