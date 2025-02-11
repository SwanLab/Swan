classdef DensityFromVademecum < VariableFromVademecum
    
   properties (Access = protected)
       fieldName = 'volume'; 
   end    
    
    methods (Access = public)
                
        function [rho,drho] = compute(obj,x)
            obj.computeParamsInfo(x);    
            obj.setValuesToInterpolator(x);
            [rho,drho] = obj.computeValues();            
        end
        
    end
    
    methods (Access = private)
        
        function [rho,drho] = computeValues(obj)
            d = squeeze(obj.values);
            [r,dr] = obj.interpolator.interpolate(d);
            rho(:,1) = r;
            drho     = dr;            
        end
        
    end
    
    
end
