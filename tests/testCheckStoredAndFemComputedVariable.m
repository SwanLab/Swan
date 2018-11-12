classdef testCheckStoredAndFemComputedVariable <  test  ...
         & testFemComputation ...
         & testLoadStoredVariable
         
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            d = numel(obj.variablesToStore);
            err = ones(d,1);
            for ivar = 1:d
                sV = obj.storedVar{ivar};
                cV = obj.computedVar{ivar};
                err(ivar) = norm(sV - cV);
            end
            hasPassed = norm(err) < 1e-6;
        end

    end
    
end

