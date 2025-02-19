classdef TensorVoigtIndexTransformer < handle

    properties (Access = protected)
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
        
        function  [istre,jstre] = voigt2tensor(obj,ind)
            [istre,jstre] = transformVoigt2Tensor(obj,ind);
        end
        
        
        function ind = tensor2Voigt(obj,istre,jstre)
            ind = obj.Tensor2Voigt(istre,jstre);
        end
        
    end
    
    methods (Access = protected, Abstract)
        loadVoigt2Tensor(obj)
        loadTensor2Voigt(obj)
    end
    
end

