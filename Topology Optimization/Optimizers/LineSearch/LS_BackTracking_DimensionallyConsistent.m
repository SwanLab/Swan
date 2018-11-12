classdef LS_BackTracking_DimensionallyConsistent < LineSearch
    methods
        function initKappa(obj,~,~,~)
            obj.kappa = 1;
        end
        
        function computeKappa(obj)
            obj.kappa = obj.kappa/obj.kfrac;
        end
    end
end

