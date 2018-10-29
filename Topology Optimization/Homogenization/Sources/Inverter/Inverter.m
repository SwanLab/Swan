classdef Inverter < handle
    
    properties (Abstract, Access = protected)
        invertedTensor
    end
    
    methods (Access = public,Static)
        
        function invertedTensor = invert(Tensor)
            factory  = InverterFactory();
            inverter = factory.create(Tensor);
            invertedTensor = inverter.getInvertedTensor();
        end
    end
    
    methods (Access = private)
        
        function T = getInvertedTensor(obj)
            T = obj.invertedTensor;
        end
        
    end
    
end

