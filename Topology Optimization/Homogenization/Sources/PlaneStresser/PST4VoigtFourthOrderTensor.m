classdef PST4VoigtFourthOrderTensor < PlaneStressTransformer
    
    
    methods (Access = protected)        
        
        function createPlaneStressTensor(obj)
            obj.psTensor = StiffnessPlaneStressVoigtTensor();
        end
        
    end
   
    
end