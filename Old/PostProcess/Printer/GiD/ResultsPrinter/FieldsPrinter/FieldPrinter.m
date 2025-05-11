classdef FieldPrinter < handle
    
    properties (Access = protected)
        fileID
        fieldValues
        fieldName
        iter
        fieldPosition
        simulationStr
        fieldRepresenter
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
        
        function printFieldLines(obj)
            obj.fieldRepresenter.printFieldLines()
        end
        
        function printResultsLineHeader(obj)
           obj.fieldRepresenter.printResultsLineHeader()
        end
        
    end
    
 end

