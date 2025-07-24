classdef FourthOrder3DTensor < AbstractTensor ...
                                        & FourthOrderDescriptor ...
                                        & TensorRepresentation ...
                                        & Elasticity3dDescriptor
    
    methods (Access = public)
        
        function obj = FourthOrder3DTensor()
        end
        
        function createRandomTensor(obj)
            obj.createRandomTensor@AbstractTensor();
        end
        
    end
    
    methods (Access = protected)
        
        function loadTensorSize(obj)
            obj.tensorSize = [3,3,3,3];
        end
    end
    
end

