classdef PlaneStressTransformer < handle
    
    
    properties (Access = protected)
        psTensor
    end
    
   
    
    methods (Access = public,Static)
        
        function psTensor = transform(Tensor)
            factory       = PlaneStressTransformerFactory();
            planeStresser = factory.create(Tensor);
            psTensor      = planeStresser.getTensor();
        end
        
    end
    
    methods (Access = private)
        
        function t = getTensor(obj)
            t = obj.psTensor;
        end
        
    end
    
    methods (Access = protected,Abstract)
        createPlaneStressTensor(obj)
    end
    
end

