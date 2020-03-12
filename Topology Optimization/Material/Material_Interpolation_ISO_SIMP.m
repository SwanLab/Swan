classdef Material_Interpolation_ISO_SIMP < Material_Interpolation
    
    properties (Access = protected)
        pExp
    end
    
    methods (Access = protected)
          
        function [k,dk] = computeKappaSymbolicFunctionAndDerivative(obj)
            [k0,k1] = obj.computeKappaLimits();            
            k = obj.interpolate(k0,k1);            
            dk = diff(k);
        end
        
        function [mu,dmu] = computeMuSymbolicFunctionAndDerivative(obj)
            [mu0,mu1] = obj.computeMuLimits();            
            mu = obj.interpolate(mu0,mu1);            
            dmu = diff(mu);
        end
        
        function f = interpolate(obj,f0,f1)
            p = obj.pExp;
            [drho0,drho1] = obj.computeDensities();            
            f = (drho0^p)*f0 + (drho1^p)*f1;
        end

        function [drho0,drho1] = computeDensities(obj)
            rho1 = obj.rho1;
            rho0 = obj.rho0;           
            rho = sym('rho','real');            
            inc  = rho1 - rho0;
            drho0 = (rho1 - rho)/inc;
            drho1 = (rho  - rho0)/inc;             
        end
               
    end
end
