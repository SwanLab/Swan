classdef LS_BackTracking_PG < LineSearch
    
    properties
        scalar_product
        kappaMultiplier
    end
    

    methods (Access = public)
        
        function obj = LS_BackTracking_PG(cParams)
            obj.scalar_product = cParams.scalarProduct;
            obj.kappaMultiplier = cParams.kappaMultiplier;
            obj.kfrac = cParams.kfrac;
            obj.kappa_min = 1e-15;
        end
        
    end 
  
end

