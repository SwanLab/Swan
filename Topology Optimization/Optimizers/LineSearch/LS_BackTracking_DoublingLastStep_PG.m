classdef LS_BackTracking_DoublingLastStep_PG < LS_BackTracking_DoublingLastStep & LS_BackTracking_PG
    methods
        function obj = LS_BackTracking_DoublingLastStep_PG(cParams)
            obj@LS_BackTracking_PG(cParams);
        end
        
        function initKappa(obj,x,g)
            if isempty(obj.kappa)
                norm_gamma = sqrt(obj.scalar_product.computeSP(x,x));
                norm_g = sqrt(obj.scalar_product.computeSP(g,g));
                obj.kappa = norm_gamma/norm_g;
            else
                obj.kappa = obj.kappaMultiplier*obj.kappa*obj.kfrac;
            end
        end
    end
end

