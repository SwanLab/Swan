classdef VoigtTensorInverter < Inverter
    
    
    methods (Access = public)
        
        function obj = VoigtTensorInverter(tensor)
            obj.compute(tensor);
        end
    end
    
    methods (Access = protected)
        
        function computeInverse(obj)
            A = obj.tensor.getValue();
            invA = inv(A);
            obj.invertedTensor.setValue(invA);
        end
        
    end
    
end

