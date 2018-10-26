classdef SymmetricTensorInverter < Inverter
    
    properties (Access = protected)
        invertedTensor
    end
    
    properties (Access = private)
       tensor 
       voigtTensor
       invVoigtTensor
    end
    
    methods (Access = public)
        
        function obj = SymmetricTensorInverter(tensor)
            obj.init(tensor)
            obj.transformTensor2Voigt()
            obj.makeInverseOfMatrix()
            obj.transformVoigt2Tensor()
        end
    end
    
    methods (Access = private)

        function init(obj,tensor)
            obj.tensor = tensor;            
        end
        
        function transformTensor2Voigt(obj)
            t = FourthOrderTensor();
            t.setValue(obj.tensor)
            tv = Tensor2VoigtConverter.convert(t); 
            obj.voigtTensor = double(tv);
        end
        
        function makeInverseOfMatrix(obj)
            obj.invVoigtTensor = inv(obj.voigtTensor);            
        end
        
        function transformVoigt2Tensor(obj)
            invT = VoigtTensor();
            invT.setValue(obj.invVoigtTensor);
            obj.invertedTensor = Voigt2TensorConverter.convert(invT); 
        end
        
        
    end
    
end

