classdef NodalFieldPrinter < handle
    
    properties (Access = protected, Abstract)
        fileID
        fieldName
        istep
        fieldType
        fieldPosition
    end
    
    methods (Access = protected)
        
        
        function printResultsLineHeader(obj)
            iD = obj.fileID;
            fN = obj.fieldName;
            is = obj.istep;
            fT = obj.fieldType;
            rL = obj.fieldPosition;
            fprintf(iD,'\nResult "%s" "%s" %.0f %s %s\n',fN,fN,is,fT,rL);
        end
    end
    
end

