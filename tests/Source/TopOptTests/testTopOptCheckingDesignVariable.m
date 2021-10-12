classdef testTopOptCheckingDesignVariable < testCheckStoredWithTopOptComputedVariable
    
  properties (Access = protected)
        variablesToStore = {'x'};
    end
       
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.computation.designVariable.value;            
        end
        
    end
end

