classdef FemMatOoHeader < handle
    
    
    methods (Access = public)
        
        function print(obj,iD)
            obj.printHeader(iD)
        end
        
    end
    
    methods (Access = private, Static)
        
        function printHeader(iD)
            fprintf(iD,'####################################################\n');
            fprintf(iD,'################# FEM-MAT-OO v.1.0 #################\n');
            fprintf(iD,'####################################################\n');
            fprintf(iD,'\n');
        end
        
    end
    
end