classdef Tensor2VoigtConverterFor3DTensors < Tensor2VoigtConverter
    
    properties (Access = protected)
        voigtTensor
        dim
        dimVoigt
        tensor
        indexTransformer
        voigtTensorSize
    end
    
    methods (Access = protected)
        
        function computeConversion(obj,tensor)
            obj.init(tensor)
            obj.createTensorVoigt();
            obj.representTensorInVoigt(obj.tensor)
        end
        
        function init(obj,tensor)
            obj.tensor = tensor;
            obj.indexTransformer = TensorVoigtIndexTransformer();
            obj.dim = size(obj.tensor,1);
            obj.dimVoigt = 6;
        end
        
        function t = createTensorVoigt(obj)
            obj.obtainVoigtTensorSize()
            isNotSymbolic = isUnit(obj.tensor);
            if isNotSymbolic
                t = zeros(obj.voigtTensorSize);
            else
                t = sym(zeros(obj.voigtTensorSize));
            end
            obj.voigtTensor = t;
        end
        
    end
    
    methods (Abstract, Access = protected)
       obtainVoigtTensorSize(obj) 
       representTensorInVoigt(obj)
    end
    
end

