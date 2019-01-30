classdef FileComparetor < handle
        
   methods (Access = public, Static)
       
       function theyAre = areFilesDifferent(fileA,fileB)           
           javaFileA = javaObject('java.io.File', fileA);
           javaFileB = javaObject('java.io.File', fileB);
           is_equal  = javaMethod('contentEquals','org.apache.commons.io.FileUtils',...
                      javaFileA, javaFileB);
           theyAre = ~is_equal;         
       end       
       
   end
    
    
end