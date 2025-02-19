classdef TensorVoigtIndexTransformer3D < TensorVoigtIndexTransformer
    
    methods (Access = protected)
        
        function loadVoigt2Tensor(obj)
            obj.Voigt2Tensor =  [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
        end
        
        function loadTensor2Voigt(obj)
            obj.Tensor2Voigt =  [1 6 5;
                                 6 2 4;
                                 5 4 3];
        end
        
    end
    
    
end