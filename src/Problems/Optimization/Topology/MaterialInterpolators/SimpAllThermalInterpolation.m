classdef SimpAllThermalInterpolation < handle
    
    properties (Access = public)
        fun
        dfun
    end

    properties (Access = private)
       f0
       f1
       ndim
    end
    
    methods (Access = public)
        
        function obj = SimpAllThermalInterpolation(cParams)
            obj.init(cParams)
            obj.computeConductivityFunctionAndDerivative();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.f0 = cParams.f0;
            obj.f1 = cParams.f1;
            switch cParams.dim
                case '2D'
                  obj.ndim = 2;
                case '3D'
                  obj.ndim = 3;
            end
        end
        
        function  computeConductivityFunctionAndDerivative(obj)
            f0V = obj.f0;
            f1V = obj.f1;
            gamma = @(rho) (1-rho)*f0V + rho*f1V;
            if obj.ndim == 2
                den = f0V + f1V;
                obj.fun  = @(rho) (gamma(rho).^2 + f0V*f1V)/den;
                obj.dfun = @(rho) (2*gamma(rho)*(f1V-f0V))/den;
            elseif obj.ndim == 3
                a1 = 3*(f1V+f0V);
                a2 = 2*f1V*f0V*(f1V+f0V);
                den = (2*f1V + f0V)*(f1V + 2*f0V);
                obj.fun  = @(rho) (-gamma(rho).^3 + a1*gamma(rho).^2 + a2)/den;
                obj.dfun = @(rho) (f1V-f0V)*(2*a1*gamma(rho) - 3*gamma(rho).^2)/den;        
            end
        end
        
    end
    
end