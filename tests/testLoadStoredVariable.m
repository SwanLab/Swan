classdef testLoadStoredVariable < handle
    
    properties (Abstract, Access = protected)
        testName
        variablesToStore
    end
    
    properties (Access = protected)
       storedVar
    end
    
    methods (Access = protected)
        
        function obj = testLoadStoredVariable()
           obj.loadStoredVariable()
        end
    end
    
    methods (Access = protected)
        function loadStoredVariable(obj)
            load_file = strcat('./tests/',obj.testName);
            load(load_file);
            for icell = 1:numel(obj.variablesToStore)
              obj.storedVar{icell} = eval(obj.variablesToStore{icell});
            end
        end
    end
    

end

