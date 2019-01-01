classdef GaussFieldPrinter < handle
    
    properties (Access = protected, Abstract)
        fileID
        fieldName
        istep
        fieldType
        fieldPosition
        simulationCase
        gaussDescriptor
    end
    
    methods (Access = protected)
        
        function printResultsLineHeader(obj)
            iD = obj.fileID;
            fN = obj.fieldName;
            sC = obj.simulationCase;
            is = obj.istep;
            fT = obj.fieldType;
            rL = obj.fieldPosition;
            gD = obj.gaussDescriptor;
            fprintf(iD,'\nResult "%s" "%s" %.0f %s %s "%s"\n',fN,sC,is,fT,rL,gD);
        end
      
    end
    
end

