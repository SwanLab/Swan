classdef DensityFromVademecum < VariableFromVademecum
    
   properties (Access = protected)
       fieldName = 'volume'; 
   end    
    
    methods (Access = public)
                
        function [rho,drho] = computeDensity(obj,x)
            obj.computeParamsInfo(x);    
            obj.setValuesToInterpolator(x);                        
            
            rho  = zeros(obj.np,1);
            drho = zeros(2*obj.np,1);
            d = squeeze(obj.values);
            [r,dr] = obj.interpolator.interpolate(d);
            rho(:,1) = r;
            drho(obj.indexMx) = dr(:,1);
            drho(obj.indexMy) = dr(:,2);
        end
        
    end
    
    
end
