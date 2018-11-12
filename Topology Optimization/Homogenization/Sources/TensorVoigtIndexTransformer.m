classdef TensorVoigtIndexTransformer < handle

    properties (Access = private)
        Voigt2Tensor
        Tensor2Voigt
    end
    
    methods (Access = public)
        
        function obj = TensorVoigtIndexTransformer()
            obj.loadVoigt2Tensor()
            obj.loadTensor2Voigt()
        end
        
        function  [istre,jstre] = transformVoigt2Tensor(obj,ind)
            istre = obj.Voigt2Tensor(ind,1);
            jstre = obj.Voigt2Tensor(ind,2);
        end
        
        function ind = transformTensor2Voigt(obj,istre,jstre)
            ind = obj.Tensor2Voigt(istre,jstre);
        end
        
    end
    
    methods (Access = private)
        
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

