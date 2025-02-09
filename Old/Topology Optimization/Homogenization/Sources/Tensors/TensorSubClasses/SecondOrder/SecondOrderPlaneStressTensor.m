classdef SecondOrderPlaneStressTensor < AbstractTensor ...
                                        & SecondOrderDescriptor ...
                                        & TensorRepresentation ...
                                        & ElasticityPlaneStressDescriptor

    methods (Access = public)

        function obj = SecondOrderPlaneStressTensor()
        end
    end

    methods (Access = protected)

        function loadTensorSize(obj)
            obj.tensorSize = [2,2];
        end
    end


end

