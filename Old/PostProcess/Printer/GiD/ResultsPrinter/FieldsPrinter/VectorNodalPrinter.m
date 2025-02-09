classdef VectorNodalPrinter < VectorPrinter
    
    methods (Access = public)
        
        function obj = VectorNodalPrinter(d)
            obj.nComp = size(d.fieldValues,2);
            obj.init(d);
            obj.createFormatString();
            obj.createFieldRepresenter();
            obj.print();
        end
    end
    
    methods (Access = private)
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                ifieldName = fieldsNames{ifield};
                fieldValue = d.(ifieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
    
        function createFieldRepresenter(obj)
            s.fileID        = obj.fileID;
            s.fieldName     = obj.fieldName;
            s.simulationStr = obj.simulationStr;
            s.iter          = obj.iter;
            s.fieldType     = obj.fieldType;
            s.fieldPosition = obj.fieldPosition;
            s.fieldValues   = obj.fieldValues;
            s.nComp         = obj.nComp;
            obj.fieldRepresenter = NodalFieldPrinter(s);
        end
        
    end
    
end

