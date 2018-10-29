classdef PlaneStressTransformer < handle
    
    
    properties (Abstract, Access = protected)
        TensorInPlaneStress
    end
    
    methods (Access = public,Static)
        
        function PlaneStressTensor = transform(Tensor)
            Factory           = PlaneStressTransformerFactory();
            PlaneStresser     = Factory.create(Tensor);
            PlaneStressTensor = PlaneStresser.getTensor();
        end
        
    end
    
    methods (Access = private)
        
        function TPS = getTensor(obj)
            TPS = double(obj.TensorInPlaneStress);
        end
        
    end
    
end

