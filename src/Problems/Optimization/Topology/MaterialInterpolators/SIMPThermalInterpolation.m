classdef SIMPThermalInterpolation < handle
    
    properties (Access = public)
        fun
        dfun
    end

    properties (Access = private)
       f0
       f1
       pExp
    end
    
    methods (Access = public)
        
        function obj = SIMPThermalInterpolation(cParams)
            obj.init(cParams)
            obj.computeConductivityFunctionAndDerivative();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.f0 = cParams.f0;
            obj.f1 = cParams.f1;
            obj.pExp = cParams.pExp;
        end
        
        function  computeConductivityFunctionAndDerivative(obj)
            f0V = obj.f0;
            f1V = obj.f1;
            p  = obj.pExp;
            obj.fun  = @(rho) (1-rho.^p)*f0V + (rho.^p)*f1V;
            obj.dfun = @(rho) -p*rho.^(p-1)*f0V + p*(rho.^(p-1))*f1V;            
        end
        
    end
    
end