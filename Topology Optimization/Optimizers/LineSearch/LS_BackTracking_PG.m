classdef LS_BackTracking_PG < LineSearch
    properties
        scalar_product
        kappaMultiplier
    end
    
    methods
        function obj = LS_BackTracking_PG(settings,epsilon)
            obj.scalar_product = ScalarProduct(settings.filename,epsilon);
            obj.kappaMultiplier = settings.kappaMultiplier;
            obj.kfrac = 2;
            obj.kappa_min = 1e-15;
        end
    end    
end

