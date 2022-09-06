classdef Tester < handle
    
    properties (Access = public, Abstract)
        testName
        calcValues
        corrValues
    end
    
    methods (Access = public)
 
        function bool = verify(obj)
            bool = 1;
            iVar = 1;
            while iVar <= size(obj.calcValues,2)
                if obj.calcValues(iVar).Matrix == obj.corrValues(iVar).Matrix
                else
                    bool = 0;
                end
                iVar = iVar + 1;
            end

        end
        
    end

    methods (Access = public, Static)
        
        function obj = create(type,initialData)
            obj = TesterFactory.create(type,initialData);            
        end
        
    end
    
end