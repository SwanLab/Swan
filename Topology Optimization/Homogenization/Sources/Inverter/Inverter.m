classdef Inverter < handle
    
    properties (Access = protected)
        invertedTensor
        tensor
    end
    
    methods (Access = public,Static)
        
        function invertedTensor = invert(Tensor)
            factory  = InverterFactory();
            inverter = factory.create(Tensor);
            invertedTensor = inverter.getInvertedTensor();
        end
    end
    
    methods (Access = protected)
        
        function compute(obj,tensor)
            obj.init(tensor)
            obj.createInvertedTensor()
            obj.computeInverse()
        end
        
    end
    
    methods (Access = private)

        function createInvertedTensor(obj)
            obj.invertedTensor = obj.tensor.clone();
        end
        
        function T = getInvertedTensor(obj)
            T = obj.invertedTensor;
        end
        
        function init(obj,t)
            obj.tensor = t;
        end
        
    end
    
    methods (Abstract, Access = protected)
       computeInverse(obj) 
    end
    
    

end

