classdef LS_BackTracking_DimensionallyConsistent_PG < LS_BackTracking_DimensionallyConsistent & LS_BackTracking_PG
    methods
        function obj = LS_BackTracking_DimensionallyConsistent_PG(settings,epsilon)
            obj@LS_BackTracking_PG(settings,epsilon);
        end
        
        function initKappa(obj,x,g)
            xNorm = obj.scalar_product.computeSP(x,x);
            gNorm = obj.scalar_product.computeSP(g,g);
            obj.kappa = xNorm/gNorm;
        end
    end
end

