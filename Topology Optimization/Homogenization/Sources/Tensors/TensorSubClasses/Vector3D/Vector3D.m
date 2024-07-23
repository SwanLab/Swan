classdef Vector3D < AbstractTensor ...
                    & FirstOrderDescriptor ...
                    & Elasticity3dDescriptor
    
    methods (Access = public)
        
        function normalize(obj)
            d = obj.getValue();
            d = d/norm(d);
            obj.setValue(d);
        end
        
    end

    methods (Access = protected)
        
        function loadOrderVariable(obj)
            obj.order = 'first';
        end
        
        function loadTensorSize(obj)
            obj.tensorSize = [3 1];
        end

    end
    
end

