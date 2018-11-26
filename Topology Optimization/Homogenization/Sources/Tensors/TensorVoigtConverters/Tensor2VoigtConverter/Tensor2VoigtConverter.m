classdef Tensor2VoigtConverter < handle
    
    properties (Abstract, Access = protected)
        voigtTensor
    end
    
    methods (Access = public,Static)
        
        function voigtTensor = convert(tensor)
            factory      = Tensor2VoigtConverterFactory();
            t2vConverter = factory.create(tensor);
            voigtTensor  = t2vConverter.getTensor();
        end
        
    end
    
    methods (Access = private)
        
        function vt = getTensor(obj)
            vt = obj.voigtTensor;
        end
        
    end
    
end

