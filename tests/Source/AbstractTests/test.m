classdef test < handle
    
    properties (Access = protected)
       FileName 
    end
    
    methods (Access = public)
        
        function checkTestPassed(obj,FileName)
            obj.FileName = FileName;
            if obj.hasPassed()                
                obj.printTestPassed()
            else
                obj.printTestNotPassed()
            end
        end
        
    end
    
    methods (Abstract, Access = protected)
        hasPassed(obj)
        printTestPassed(obj)
        printTestNotPassed(obj)
    end
    
end

