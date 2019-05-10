classdef LS_BackTracking_PG < LineSearch
    properties
        scalar_product
        kappaMultiplier
    end
    

    methods
        function obj = LS_BackTracking_PG(cParams,epsilon)
            obj.scalar_product = ScalarProduct(cParams.scalarProductSettings);
            obj.kappaMultiplier = cParams.kappaMultiplier;
            obj.kfrac = cParams.kfrac;
            obj.kappa_min = 1e-15;
        end
    end 
    
    
end

