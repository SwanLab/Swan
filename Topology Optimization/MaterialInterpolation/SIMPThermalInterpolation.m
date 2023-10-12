classdef SIMPThermalInterpolation < handle
       
    properties (Access = private)
       alpha0
       alpha1
    end
    
    properties (Access = private)
        pExp
        alpha
        dalpha
    end
    
    methods (Access = public)
        
        function obj = SIMPThermalInterpolation(cParams)
            obj.init(cParams)
            obj.computeConductivityFunctionAndDerivative();
        end
        
        function mp = computeMatProp(obj,rho)
            mp.alpha  = obj.alpha(rho);
            mp.dalpha = obj.dalpha(rho);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.alpha0 = cParams.alpha0;
            obj.alpha1 = cParams.alpha1;
            obj.pExp = 3;            
        end
        
        function  computeConductivityFunctionAndDerivative(obj)
            a0 = obj.alpha0;
            a1 = obj.alpha1;
            p  = obj.pExp;
            obj.alpha = @(rho) (1-rho.^p)*a0 + (rho.^p)*a1;
            obj.dalpha = @(rho) -p*rho.^(p-1)*a0 + p*(rho.^(p-1))*a1;            
        end
        
    end
    
end