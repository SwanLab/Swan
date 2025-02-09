classdef StrainDescriptor < FieldNameDescriptor
        
    methods (Access = protected)
        
        function loadFieldVariable(obj)
            obj.fieldName = 'strain';            
        end
        
    end
    
end

