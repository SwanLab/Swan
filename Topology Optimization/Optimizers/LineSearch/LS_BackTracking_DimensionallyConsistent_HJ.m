classdef LS_BackTracking_DimensionallyConsistent_HJ < LS_BackTracking_DimensionallyConsistent & LS_BackTracking_HJ    
    methods
        function obj = LS_BackTracking_DimensionallyConsistent_HJ(settings)
            obj@LS_BackTracking_HJ(settings);
        end
        
        function initKappa(obj,~,~,~)
            obj.kappa = 1;
            obj.HJiter = obj.HJiter0;
        end
        
        function computeKappa(obj)
            if obj.HJiter > obj.HJiter_min
                obj.HJiter = round(obj.HJiter/obj.kfrac);
            else
                obj.kappa = obj.kappa/obj.kfrac;
            end
        end
    end
end

