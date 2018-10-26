classdef Voigt2TensorConverterFor3DTensors < Voigt2TensorConverter
    
    properties (Access = protected)
        voigtTensor
        dim
        dimVoigt
        tensor
        indexTransformer
        tensorSize
    end
    
    methods (Access = protected)
        
        function computeConversion(obj,tensor)
            obj.init(tensor)
            obj.createTensor();
            obj.representVoigtInTensor(obj.voigtTensor)
        end
        
        function init(obj,tensor)
            obj.voigtTensor = tensor;
            obj.indexTransformer = TensorVoigtIndexTransformer();
            obj.dim = 3;
            obj.dimVoigt = 6;
        end
        
        function t = createTensor(obj)
            obj.obtainTensorSize()
            isNotSymbolic = isUnit(obj.voigtTensor);
            if isNotSymbolic
                t = zeros(obj.tensorSize);
            else
                t = sym(zeros(obj.tensorSize));
            end
            obj.tensor = t;
        end
        
    end
    
    methods (Abstract, Access = protected)
       obtainTensorSize(obj) 
       representVoigtInTensor(obj)
    end
    
end
