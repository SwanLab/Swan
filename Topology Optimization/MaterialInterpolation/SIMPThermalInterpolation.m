classdef SIMPThermalInterpolation < L2Function
    
    properties (Access = public)
        ndimf
    end

    properties (Access = private)
       alpha0
       alpha1
       density       
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
        
        function alpha = evaluate(obj,xV)
            rho = obj.density.evaluate(xV);
            alpha  = obj.alpha(rho);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.alpha0 = cParams.alpha0;
            obj.alpha1 = cParams.alpha1;
            obj.density = cParams.density; 
            obj.pExp = 3;        
            obj.ndimf = 1;
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