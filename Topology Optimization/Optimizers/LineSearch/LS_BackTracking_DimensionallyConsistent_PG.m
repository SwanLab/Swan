classdef LS_BackTracking_DimensionallyConsistent_PG < LS_BackTracking_DimensionallyConsistent & LS_BackTracking_PG
    methods
        function obj = LS_BackTracking_DimensionallyConsistent_PG(settings,epsilon)
            obj@LS_BackTracking_PG(settings,epsilon);
        end
        
        function initKappa(obj,x,gradient)
            norm_gamma = sqrt(obj.scalar_product.computeSP(x,x));
            norm_g = sqrt(obj.scalar_product.computeSP(gradient,gradient));
            obj.kappa = norm_gamma/norm_g;
        end
    end
end

