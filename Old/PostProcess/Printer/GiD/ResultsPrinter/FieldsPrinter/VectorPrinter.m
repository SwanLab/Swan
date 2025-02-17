classdef VectorPrinter < FieldPrinter
    
    properties (Access = protected)
        fieldType = 'Vector';
        fieldComponentName
        formatString
        nComp
    end
    
    methods (Access = protected)
        
        function print(obj)
            obj.printResultsLineHeader();
            obj.printComponentNamesLine();
            obj.printValuesLine();
            obj.printFieldLines();
            obj.printEndValuesLine();
        end
        
        function printComponentNamesLine(obj)
            iD = obj.fileID;
            fC = obj.fieldComponentName;
            d = obj.nComp;
            fcV = cell(1,d);
            [fcV{:}] = deal(fC);
            fprintf(iD,obj.formatString,fcV{:});
        end
        
        function createFormatString(obj)
            switch obj.nComp
                case 2
                    obj.formatString = 'ComponentNames  "%sx", "%sy"\n';
                case 3
                    obj.formatString = 'ComponentNames  "%sx", "%sy", "%sz"\n';
                case 4
                    obj.formatString = 'ComponentNames  "%sx", "%sy", "%sxy", "%sz"\n';
                case 6
                    obj.formatString = 'ComponentNames  "%sx", "%sy", "%sz", "%sxy", "%syz", "%sxz"\n';
            end
        end
        

    end
    
end

