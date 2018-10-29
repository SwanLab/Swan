classdef VoigtTensorInverter < Inverter
    
    properties (Access = protected)
        invertedTensor
    end
    
    methods (Access = public)
        
        function obj = VoigtTensorInverter(Tensor)
            obj.computeInverse(Tensor)
        end
    end
    
    methods (Access = private)
       
        function computeInverse(obj,Tensor)
            obj.invertedTensor = inv(Tensor);
        end
        
    end
    
end

