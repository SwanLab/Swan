classdef LS_BackTracking_DoublingLastStep_HJ < LS_BackTracking_DoublingLastStep & LS_BackTracking_HJ
    methods
        function obj = LS_BackTracking_DoublingLastStep_HJ(settings)
            obj@LS_BackTracking_HJ(settings);
        end
        
        function initKappa(obj,~,~,~)
            if obj.HJiter > obj.HJiter_min
                obj.HJiter = obj.HJiter*obj.kfrac;
            else
                if obj.kappa < obj.kappa_min
                    obj.kappa = obj.kappa_min*obj.kfrac;
                end
                obj.kappa = obj.kappa*obj.kfrac;
            end
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

