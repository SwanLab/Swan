classdef FieldPrinter < handle
    
    properties (Access = protected)
        fileID
        fieldValues
        fieldName
        iter
        fieldPosition
        simulationStr
    end

    properties (Access = protected, Abstract)
       fieldType
    end
    
    methods (Access = protected)
        
        function printValuesLine(obj)
            iD = obj.fileID;
            fprintf(iD,'Values\n');
        end
        
        function printEndValuesLine(obj)
            iD = obj.fileID;
            fprintf(iD,'End Values\n');
        end
    end
    
    methods (Access = protected, Abstract)
        printFieldLines(obj)
    end
    
end

