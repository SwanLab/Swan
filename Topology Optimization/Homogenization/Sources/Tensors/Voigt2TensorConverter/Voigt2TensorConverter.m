classdef Voigt2TensorConverter < handle
    
    properties (Abstract, Access = protected)
        tensor
    end
    
  methods (Access = public,Static)
        
        function tensor = convert(voigtTensor)
            factory      = Voigt2TensorConverterFactory();
            v2tConverter = factory.create(voigtTensor);
            tensor       = v2tConverter.getTensor();
        end
       
  end
            
    methods (Access = private)
        
        function t = getTensor(obj)
            t = obj.tensor;
        end
        
    end
    
end

