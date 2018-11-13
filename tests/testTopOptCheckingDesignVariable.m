classdef testTopOptCheckingDesignVariable < testCheckStoredWithTopOptComputedVariable
    
  properties (Access = protected)
        computedVar
        variablesToStore = {'x'};
    end
       
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.topOpt.x;            
        end
        
    end
end

