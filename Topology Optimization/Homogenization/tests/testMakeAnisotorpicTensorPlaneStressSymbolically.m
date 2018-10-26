classdef testMakeAnisotorpicTensorPlaneStressSymbolically < test
    
    properties (Access = protected)
        ToCheckTensor
        CheckedTensor
    end
    
    methods
        
        function obj = testMakeAnisotorpicTensorPlaneStressSymbolically()
           Tensor = FourthOrderTensor();
           Tensor.createRandomTensor();
           Tensor.computeTensorVoigt();
           
           TensorTransformer = PlaneStressVoigtTensorSymbolicallyTransformer(Tensor.tensorVoigt);
           
           obj.ToCheckTensor = TensorTransformer.getValue();           
           obj.CheckedTensor = PlaneStressTransformer.transform(Tensor.tensorVoigt);
        end
        
    end
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            hasPassed = norm(double(obj.ToCheckTensor(:)) - double(obj.CheckedTensor(:))) < 1e-13;
        end
    end
end


