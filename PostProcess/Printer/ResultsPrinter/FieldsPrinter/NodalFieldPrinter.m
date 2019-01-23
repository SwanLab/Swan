classdef NodalFieldPrinter < handle
    
    properties (Access = protected, Abstract)
        fileID
        fieldName
        iter
        fieldType
        fieldPosition
        simulationStr
    end
    
    methods (Access = protected)
        
        
        function printResultsLineHeader(obj)
            iD = obj.fileID;
            fN = obj.fieldName;
            sC = obj.simulationStr;
            it = obj.iter;
            fT = obj.fieldType;
            rL = obj.fieldPosition;
            fprintf(iD,'\nResult "%s" "%s" %.0f %s %s\n',fN,sC,it,fT,rL);
        end
    end
    
end

