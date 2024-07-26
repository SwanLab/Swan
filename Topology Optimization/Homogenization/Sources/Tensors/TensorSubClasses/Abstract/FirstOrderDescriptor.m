classdef FirstOrderDescriptor < OrderDescriptor
    

    methods (Access = protected)
        
        function loadOrderVariable(obj)
            obj.order = 'first';
        end
        
    end
    
end

