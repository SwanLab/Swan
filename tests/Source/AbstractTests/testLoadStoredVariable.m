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
    
    methods (Access = private)
        function loadStoredVariable(obj)
            file2load = obj.testName;
            load(file2load);
            for icell = 1:numel(obj.variablesToStore)
              obj.storedVar{icell} = eval(obj.variablesToStore{icell});
            end
        end
    end
end

