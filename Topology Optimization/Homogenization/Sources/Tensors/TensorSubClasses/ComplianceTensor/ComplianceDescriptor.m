classdef ComplianceDescriptor < FieldNameDescriptor
        
    methods (Access = protected)
        
        function loadFieldVariable(obj)
            obj.fieldName = 'compliance';
        end
        
    end
    
end

