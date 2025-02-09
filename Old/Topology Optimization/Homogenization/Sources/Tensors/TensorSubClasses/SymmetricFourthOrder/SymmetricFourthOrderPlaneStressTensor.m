classdef SymmetricFourthOrderPlaneStressTensor <  AbstractTensor ...
                                              & FourthOrderDescriptor ...
                                              & TensorRepresentation ...
                                              & ElasticityPlaneStressDescriptor

    methods (Access = public)
        
        function obj = SymmetricFourthOrderPlaneStressTensor()
        end
        
        function d = getDimension(obj)
            d = obj.dim;
        end
        
    end
    
    methods (Access = protected)
        
        function loadTensorSize(obj)
            obj.tensorSize = [2,2,2,2];
        end
    end
    
end

