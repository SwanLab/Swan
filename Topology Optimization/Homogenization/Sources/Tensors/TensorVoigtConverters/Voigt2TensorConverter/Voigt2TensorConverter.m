classdef Voigt2TensorConverter < handle
    
    properties (Access = protected)
        voigtTensor
        tensor
        tensorSize
    end
    
    methods (Access = public,Static)
        
        function tensor = convert(voigtTensor)
            factory      = Voigt2TensorConverterFactory();
            v2tConverter = factory.create(voigtTensor);
            tensor       = v2tConverter.getTensor();
        end
        
    end
    
    methods (Access = protected)
        function computeConversion(obj,tensor)
            obj.init(tensor)
            obj.createTensor();
            obj.representVoigtInTensor()
        end
  
    end
    
    methods (Access = private)
        
        function t = getTensor(obj)
            t = obj.tensor;
        end
         
        function init(obj,tensor)
            obj.voigtTensor = tensor;
            obj.createIndexTransformer();
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
       createIndexTransformer(obj)
       selectTensorClass(obj)
    end
    
end

