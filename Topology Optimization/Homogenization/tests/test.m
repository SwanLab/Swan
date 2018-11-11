classdef test < handle

    properties
    end

    methods
  
        function checkTestPassed(obj,FileName)
            if obj.hasPassed()
                cprintf('green',strcat(FileName,' PASSED\n'));
            else
                cprintf('err',strcat(FileName,' FAILED\n'));
            end
        end
        

    end
    
    methods (Abstract, Access = protected)
        hasPassed(obj)        
    end


end

