classdef LS_BackTracking_DoublingLastStep_PG < LS_BackTracking_DoublingLastStep & LS_BackTracking_PG
    methods
        function obj = LS_BackTracking_DoublingLastStep_PG(settings,epsilon)
            obj@LS_BackTracking_PG(settings,epsilon);
        end
        
        function initKappa(obj,x,gradient)
            if isempty(obj.kappa)
                norm_gamma = sqrt(obj.scalar_product.computeSP(x,x));
                norm_g = sqrt(obj.scalar_product.computeSP(gradient,gradient));
                obj.kappa = norm_gamma/norm_g;
            else
                obj.kappa = obj.kappaMultiplier*obj.kappa*obj.kfrac;
            end
        end
    end
end

