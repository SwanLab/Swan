classdef ScalarNodalPrinter < ScalarPrinter
    
    properties (Access = private)
        nDecimalPositions = 12
    end
    
    methods (Access = public)
        
        function obj = ScalarNodalPrinter(d)
            obj.init(d);
            obj.createFieldRepresenter();                        
            obj.print();
        end
        
    end
        
    methods (Access = private)
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                fieldName = fieldsNames{ifield};
                fieldValue = d.(fieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end        
        
        function roundFieldValues(obj)
            obj.fieldValues = round(obj.fieldValues,obj.nDecimalPositions);
        end
        
        function createFieldRepresenter(obj)
            s.fileID        = obj.fileID;
            s.fieldName     = obj.fieldName;
            s.simulationStr = obj.simulationStr;
            s.iter          = obj.iter;
            s.fieldType     = obj.fieldType;
            s.fieldPosition = obj.fieldPosition;
            s.fieldValues   = obj.fieldValues;
            s.nComp         = 1;            
            obj.fieldRepresenter = NodalFieldPrinter(s);
        end                
        
    end
    
end

