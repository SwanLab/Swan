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
            fprintf(iD,'####################################################\n');
            fprintf(iD,'################# FEM-MAT-OO v.1.0 #################\n');
            fprintf(iD,'####################################################\n');
            fprintf(iD,'\n');
        end
        
    end
    
end