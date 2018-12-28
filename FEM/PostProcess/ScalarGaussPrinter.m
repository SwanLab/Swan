classdef ScalarGaussPrinter < FieldPrinter ...
                       & GaussFieldPrinter
    
    properties (Access = protected)
        gaussDescriptor
    end
                   
    methods (Access = public)
        
        function obj = ScalarGaussPrinter(d)
            obj.fieldType = 'Scalar';
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
            fV = obj.fieldValues;
            for inode = 1:length(fV)
                fprintf(iD,'%6.0f ',inode);
                fprintf(iD,'%12.5d ',fV(inode));
                fprintf(iD,'\n');
            end
        end
        
    end
    
end

