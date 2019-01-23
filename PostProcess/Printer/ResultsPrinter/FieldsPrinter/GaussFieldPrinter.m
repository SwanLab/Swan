classdef GaussFieldPrinter < handle
    
    properties (Access = protected, Abstract)
        fileID
        fieldName
        iter
        fieldType
        fieldPosition
        simulationStr
        gaussDescriptor
    end
    
    methods (Access = protected)
        
        function printResultsLineHeader(obj)
            iD = obj.fileID;
            fN = obj.fieldName;
            sC = obj.simulationStr;
            is = obj.iter;
            fT = obj.fieldType;
            rL = obj.fieldPosition;
            gD = obj.gaussDescriptor;
            fprintf(iD,'\nResult "%s" "%s" %.0f %s %s "%s"\n',fN,sC,is,fT,rL,gD);
        end
      
    end
    
end

