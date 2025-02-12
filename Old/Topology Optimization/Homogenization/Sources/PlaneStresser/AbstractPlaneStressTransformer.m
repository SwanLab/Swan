classdef AbstractPlaneStressTransformer < handle
    
    
    properties (Abstract, Access = protected)
        TensorInPlaneStress
    end
    
    
    methods (Access = private)
        
        function TPS = getTensor(obj)
            TPS = obj.TensorInPlaneStress;
        end
        
    end
    
    
    methods (Access = public,Static)
        
        function PlaneStressTensor = transform(Tensor)
            Factory           = PlaneStressTransformerFactor();
            PlaneStresser     = Factory.create(Tensor);
            PlaneStressTensor = PlaneStresser.getTensor();
        end
        
        
    end
    
end

