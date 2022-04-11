classdef Tester < handle
    
    properties (Access = protected, Abstract)
        testName
        calcValues
        corrValues
    end
    
    methods (Access = public, Static)
        
        function obj = create(type,initialData)
            obj = TesterFactory.create(type,initialData);            
        end

    end
    
    methods (Access = protected)
 
        function verify(obj)
            bool = 1;
            iVar = 1;
            while iVar <= size(obj.calcValues,2)
                if obj.calcValues(iVar).Matrix == obj.corrValues(iVar).Matrix
                else
                    bool = 0;
                end
                iVar = iVar + 1;
            end
            if bool == 1
                cprintf('green',['Test pass. ', obj.testName, ' working properly.\n']);
            else
                cprintf('red', ['Test NO pass. ', obj.testName, ' failed.\n']);
            end
            
        end
        
    end
    
end