classdef SecondOrderDescriptor < OrderDescriptor
    


    methods (Access = protected)
        
        function loadOrderVariable(obj)
            obj.order = 'second';            
        end
        
    end
    
end

