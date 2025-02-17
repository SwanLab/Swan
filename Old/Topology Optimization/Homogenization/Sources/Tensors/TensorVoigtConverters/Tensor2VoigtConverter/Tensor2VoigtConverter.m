classdef Tensor2VoigtConverter < handle
    
    properties (Access = protected)
        voigtTensor
        tensor
        voigtTensorSize
    end
    
    methods (Access = public,Static)
        
        function voigtTensor = convert(tensor)
            factory      = Tensor2VoigtConverterFactory();
            t2vConverter = factory.create(tensor);
            voigtTensor  = t2vConverter.getTensor();
        end
        
    end
    
    methods (Access = protected)
        
        function computeConversion(obj,tensor)
            obj.init(tensor)
            obj.createTensorVoigt();
            obj.representTensorInVoigt()
        end
        
    end
    
    
    methods (Access = private)
        
        function vt = getTensor(obj)
            vt = obj.voigtTensor;
        end
        
        function init(obj,tensor)
            obj.tensor = tensor;
            obj.createIndexTransformer();
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
        
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isSymbolic(v)
            itIs = isa(v,'sym');
        end
        
    end
    
    methods (Abstract, Access = protected)
        representTensorInVoigt(obj)
        selectVoigtTensorClass(obj)
        createIndexTransformer(obj)
    end
    
end

