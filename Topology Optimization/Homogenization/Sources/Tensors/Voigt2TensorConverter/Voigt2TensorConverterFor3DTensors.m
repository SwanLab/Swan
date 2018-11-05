classdef Voigt2TensorConverterFor3DTensors < Voigt2TensorConverter
    
    properties (Access = protected)
        voigtTensor
        tensor
        indexTransformer
        tensorSize
    end
    
    methods (Access = protected)
        
        function computeConversion(obj,tensor)
            obj.init(tensor)
            obj.createTensor();
            obj.representVoigtInTensor()
        end
        
        function init(obj,tensor)
            obj.voigtTensor = tensor;
            obj.indexTransformer = TensorVoigtIndexTransformer();
        end
        
        function createTensor(obj)
            obj.selectTensorClass();
            obj.initializeTensor(); 
        end
        
        function initializeTensor(obj)
            vt = obj.voigtTensor.getValue();
            s = obj.tensor.getTensorSize();
            t = zeros(s);
            if obj.isSymbolic(vt)
                t = sym(t);
            end
            obj.tensor.setValue(t);
        end
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isSymbolic(v)
            itIs = isa(v,'sym');
        end
        
    end
    
    methods (Abstract, Access = protected)
       representVoigtInTensor(obj)
       selectTensorClass(obj)
    end
    
end
