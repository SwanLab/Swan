classdef StressDescriptor < FieldNameDescriptor
    
    methods (Access = protected)
        
        function loadFieldVariable(obj)
            obj.fieldName = 'stress';
        end
        
    end
    
end

