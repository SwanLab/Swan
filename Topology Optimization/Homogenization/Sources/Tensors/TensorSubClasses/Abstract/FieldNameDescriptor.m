classdef FieldNameDescriptor < handle
    
    properties (Access = protected)
        fieldName
    end
    
    methods (Access = public)
        
        function obj = FieldNameDescriptor()
            obj.loadFieldVariable()
        end
        
        function o = getFieldName(obj)
            o = obj.fieldName;
        end
        
    end
    
    methods (Abstract, Access = protected)
        loadFieldVariable(obj)
    end
    
end

