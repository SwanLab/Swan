classdef testStoredComputedChecker < handle
    
    properties (Abstract, Access = protected)
        computedVar
        variablesToStore
        storedVar
        error
    end
    
    properties (Access = protected)
        tol = 1e-6;
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            d = numel(obj.variablesToStore);
            err = ones(d,1);
            for ivar = 1:d
                sV = obj.storedVar{ivar};
                cV = obj.computedVar{ivar};
                err(ivar) = norm(sV - cV)/norm(sV);
            end
            obj.error = norm(err);
        end
        
    end
    

end

