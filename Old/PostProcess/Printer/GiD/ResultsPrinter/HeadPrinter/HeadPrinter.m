classdef HeadPrinter < handle
    
    properties (Access = protected)
       fileID 
    end

    methods (Access = protected)
        
        function printInitialLine(obj)
            fprintf(obj.fileID,'GiD Post Results File 1.0\n\n');
        end
        
        function printFemMatOoHeader(obj)
            iD = obj.fileID;
            h = FemMatOoHeader();
            h.print(iD);
        end
        
    end
    
end