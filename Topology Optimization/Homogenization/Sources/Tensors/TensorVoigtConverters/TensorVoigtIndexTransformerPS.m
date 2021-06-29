classdef TensorVoigtIndexTransformerPS < TensorVoigtIndexTransformer
    
    
    
    methods (Access = protected)
        
        function loadVoigt2Tensor(obj)
            obj.Voigt2Tensor =  [1 1; 2 2; 1 2];
        end
        
        function loadTensor2Voigt(obj)
            obj.Tensor2Voigt =  [1 3;
                                 3 2];
        end
        
    end
    
    
end