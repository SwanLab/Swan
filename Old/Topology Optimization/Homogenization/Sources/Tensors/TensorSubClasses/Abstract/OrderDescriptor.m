classdef OrderDescriptor < handle
    
    properties (Access = protected)
        order
    end
        
    methods (Access = public)
        
        function obj = OrderDescriptor()
            obj.loadOrderVariable()
        end
        
        function o = getOrder(obj)
            o = obj.order;
        end
        
    end
    
    methods (Access = protected, Abstract)
        loadOrderVariable(obj)
    end
    
end

