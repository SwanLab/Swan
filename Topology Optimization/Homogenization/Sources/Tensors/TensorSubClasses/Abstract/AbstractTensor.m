classdef AbstractTensor < handle
   
    properties (Access = protected)
        tensorValue
        tensorSize
    end
        
    methods (Access = public)
        
        function obj = AbstractTensor()
            obj.loadTensorSize()
        end
       
        function T = getValue(obj)
           T = obj.tensorValue; 
        end
        
        function setValue(obj,T)
            obj.tensorValue = T;
        end
        
        function s = getTensorSize(obj)
            s = obj.tensorSize;
        end
        
        function createRandomTensor(obj)
            obj.tensorValue = rand(obj.tensorSize);
        end
        
        function t = clone(obj)
           t = eval(class(obj));
        end
        
    end
    
    methods (Abstract,Access = protected)
        loadTensorSize(obj)
    end

end

