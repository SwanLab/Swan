classdef testNotShowingError < test
    
    methods (Access = protected)
        
        function printTestPassed(obj)
           cprintf('green',obj.FileName);                        
           cprintf('green',' PASSED\n');
        end
        
        function printTestNotPassed(obj)
            cprintf('red',obj.FileName);                        
            cprintf('red',' FAILED\n');
        end
    end
end

