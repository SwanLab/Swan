classdef NodalFieldPrinter < handle
    
    properties (Access = protected, Abstract)
        fileID
        fieldName
        istep
        fieldType
        fieldPosition
        simulationCase
    end
    
    methods (Access = protected)
        
        
        function printResultsLineHeader(obj)
            iD = obj.fileID;
            fN = obj.fieldName;
            sC = obj.simulationCase;
            is = obj.istep;
            fT = obj.fieldType;
            rL = obj.fieldPosition;
            fprintf(iD,'\nResult "%s" "%s" %.0f %s %s\n',fN,sC,is,fT,rL);
        end
    end
    
end

