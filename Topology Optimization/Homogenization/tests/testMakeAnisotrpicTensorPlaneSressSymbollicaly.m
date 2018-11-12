classdef testMakeAnisotrpicTensorPlaneSressSymbollicaly < test
    
    properties (Access = protected)
        ToCheckTensor
        CheckedTensor
    end
    
    methods
        
        function obj = testMakeAnisotrpicTensorPlaneSressSymbollicaly()
           Tensor = fourthOrderTensor();
           Tensor.createRandomTensor();
           Tensor.computeTensorVoigt();
           
           TensorTransformer = PlaneStressVoigtTensorSymbolicallyTransformer(Tensor.tensorVoigt);
           
           obj.ToCheckTensor = TensorTransformer.tensorVoigtInPlaneStress;           
           obj.CheckedTensor = Tensor.transform3D_2_PlaneStressInVoigt(Tensor.tensorVoigt);
        end
        
    end
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            hasPassed = norm(double(obj.ToCheckTensor(:)) - double(obj.CheckedTensor(:))) < 1e-13;
        end
    end
end


