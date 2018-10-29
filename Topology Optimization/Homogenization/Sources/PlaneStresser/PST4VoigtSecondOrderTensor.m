classdef PST4VoigtSecondOrderTensor < PlaneStressTransformer
    
    properties (Access = protected)
       TensorInPlaneStress
    end
    
    properties (Access = private)
        Tensor
    end
    
    methods (Access = public)
        
        function obj = PST4VoigtSecondOrderTensor(Tensor)
            obj.Tensor = Tensor;
            obj.computePlaneStressTensor()
        end
        
    end
    
    methods (Access = private)

        function computePlaneStressTensor(obj)
            PSIndex = PlaneStressIndex();
            InPlaneIndex = PSIndex.getInPlaneIndex();
            obj.TensorInPlaneStress = obj.Tensor(InPlaneIndex);
        end
        
    end
        
end

