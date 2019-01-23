classdef LS_BackTracking_DoublingLastStep < LineSearch
    methods
        function initKappa(obj,~,~,~)
            if obj.kappa < obj.kappa_min
                obj.kappa = obj.kappa_min*obj.kfrac;
            end
            obj.kappa = obj.kappa*obj.kfrac;
        end
        
        function computeKappa(obj)
            obj.kappa = obj.kappa/obj.kfrac;
        end
    end
end

