classdef FileComparator < handle
    
    methods (Access = public)
        
        function theyAre = areFilesDifferent(obj,fileA,fileB)
            
            A = obj.getContent(fileA);
            B = obj.getContent(fileB);
            
            is_equal = obj.isContentEqual(A,B);
            
            theyAre = ~is_equal;
        end
        
    end
    
    methods (Access = private, Static)
        
        function a = getContent(file)
            fid = fopen(file);
            
            a = textscan(fid,'%s','delimiter','\n');
            a = a{1};
            
            fclose(fid);
        end
        
        function itIs = isContentEqual(A,B)
            itIs = all(strcmp(A,B));
        end
        
    end
    
end