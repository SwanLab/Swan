classdef StiffnessDescriptor < FieldNameDescriptor
        
    methods (Access = protected)
        
        function loadFieldVariable(obj)
            obj.fieldName = 'stiffness';            
        end
        
    end
    
end

