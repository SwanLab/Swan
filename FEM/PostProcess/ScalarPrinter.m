classdef ScalarPrinter < FieldPrinter ...
                       & NodalFieldPrinter
    
    
    methods (Access = public)
        
        function obj = ScalarPrinter(fileID,fieldValues,fieldName,istep,fieldPosition)
            obj.fieldType = 'Scalar';
            obj.init(fileID,fieldValues,fieldName,istep,fieldPosition);
            obj.print();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,fileID,fieldValues,fieldName,istep,fieldPosition)
            obj.fileID             = fileID;
            obj.fieldValues        = fieldValues;
            obj.fieldName          = fieldName;
            obj.istep              = istep;
            obj.fieldPosition      = fieldPosition;
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

