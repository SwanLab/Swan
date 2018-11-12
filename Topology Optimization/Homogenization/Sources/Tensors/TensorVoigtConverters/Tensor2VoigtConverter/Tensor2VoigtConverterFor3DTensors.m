classdef Tensor2VoigtConverterFor3DTensors < Tensor2VoigtConverter
    
    properties (Access = protected)
        voigtTensor
        tensor
        indexTransformer
        voigtTensorSize
    end
    
    methods (Access = protected)
        
        function computeConversion(obj,tensor)
            obj.init(tensor)
            obj.createTensorVoigt();
            obj.representTensorInVoigt()
        end
        
        function init(obj,tensor)
            obj.tensor = tensor;
            obj.indexTransformer = TensorVoigtIndexTransformer();
        end
        
        function createTensorVoigt(obj)
            obj.selectVoigtTensorClass();
            obj.initializeVoigtTensor();            
        end
        
        function initializeVoigtTensor(obj)
            t = obj.tensor.getValue();
            s = obj.voigtTensor.getTensorSize();
            vt = zeros(s);
            if obj.isSymbolic(t)
                vt = sym(vt);
            end
            obj.voigtTensor.setValue(vt);
        end
        
        function itIs = is3D(obj)
            itIs = strcmp(obj.tensor.getElasticityCase(),'3D');
        end
        
        function itIs = isPlaneStress(obj)
            itIs = strcmp(obj.tensor.getElasticityCase(),'planeStress');
        end
        
    end
    
    methods (Access = private, Static)

        function itIs = isSymbolic(v)
            itIs = isa(v,'sym');
        end
        
    end
       
    methods (Abstract, Access = protected)
       representTensorInVoigt(obj)
       selectVoigtTensorClass(obj)
    end
    
end

