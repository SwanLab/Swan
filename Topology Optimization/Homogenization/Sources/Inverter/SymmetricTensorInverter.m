classdef SymmetricTensorInverter < Inverter
    
    
    properties (Access = private)
        voigtTensor
        invVoigtTensor
    end
    
    methods (Access = public)
        
        function obj = SymmetricTensorInverter(tensor)
            obj.compute(tensor);
        end
    end
    
    methods (Access = protected)
        
        function computeInverse(obj)
            obj.transformTensor2Voigt()
            obj.makeInverseOfVoigtTensor()
            obj.transformVoigt2Tensor()
        end
    end
    
    methods (Access = private)
        
        
        function transformTensor2Voigt(obj)
            obj.voigtTensor = Tensor2VoigtConverter.convert(obj.tensor);
        end
        
        function makeInverseOfVoigtTensor(obj)
            obj.invVoigtTensor = Inverter.invert(obj.voigtTensor);
        end
        
        function transformVoigt2Tensor(obj)
            invT = obj.invVoigtTensor;
            obj.invertedTensor = Voigt2TensorConverter.convert(invT);
        end
        
        
    end
    
end

