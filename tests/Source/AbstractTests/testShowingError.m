classdef testShowingError < test
   
    properties (Access = protected)
       error
    end
    
    properties (Abstract, Access = protected)
        tol
    end
        
    methods (Access = protected)
        function printTestPassed(obj)
           cprintf('green',obj.FileName);                                    
           cprintf('green',' PASSED.');
           cprintf('black',['Error: ',num2str(obj.error),'\n']);
        end
        
        function printTestNotPassed(obj)
            cprintf('red',obj.FileName);                        
            cprintf('red',' FAILED.');
            cprintf('red',['Error: ',num2str(obj.error),'\n']);
        end
        
        function hasPassed = hasPassed(obj)
            obj.computeError()
            hasPassed = obj.error < obj.tol();
        end
    end
    
    methods (Abstract, Access = protected)
        computeError(obj)        
    end
end

