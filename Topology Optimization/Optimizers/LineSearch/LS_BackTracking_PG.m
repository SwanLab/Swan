classdef LS_BackTracking_PG < LineSearch
    properties
        scalar_product
        kappaMultiplier
    end
    
    methods
        function obj = LS_BackTracking_PG(settings,epsilon)
            spS.filename = settings.filename;
            spS.epsilon  = epsilon;
            obj.scalar_product = ScalarProduct(spS);
            obj.kappaMultiplier = settings.kappaMultiplier;
            obj.kfrac = 2;
            obj.kappa_min = 1e-15;
        end
    end    
end

