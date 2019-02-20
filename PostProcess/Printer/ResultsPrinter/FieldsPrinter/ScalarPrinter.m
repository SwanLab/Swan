classdef ScalarPrinter < FieldPrinter ...
        & NodalFieldPrinter
    
    properties (Access = protected)
        fieldType = 'Scalar';
    end
    
    properties (Access = private)
        nDecimalPositions = 12
    end
    
    methods (Access = public)
        
        function obj = ScalarPrinter(d)
            obj.init(d);
            obj.print();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                fieldName = fieldsNames{ifield};
                fieldValue = d.(fieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
        function print(obj)
            obj.printResultsLineHeader()
            obj.printValuesLine();
            obj.printFieldLines();
            obj.printEndValuesLine();
        end
        
        function printFieldLines(obj)
            iD = obj.fileID;
            obj.roundFieldValues();
            fV = obj.fieldValues;
            for inode = 1:length(fV)
                fprintf(iD,'%6.0f ',inode);
                fprintf(iD,'%12.5d ',fV(inode));
                fprintf(iD,'\n');
            end
        end
        
    end
    
    methods (Access = private)
        
        function roundFieldValues(obj)
            obj.fieldValues = round(obj.fieldValues,obj.nDecimalPositions);
        end
        
    end
    
end

