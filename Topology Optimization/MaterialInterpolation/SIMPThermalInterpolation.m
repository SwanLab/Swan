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
        
        function a = getFunction(obj)
            s.fHandle = cParams.fHandle;
            s.ndimf   = cParams.ndimf;
            s.mesh    = cParams.mesh;            
            a  = obj.alpha;
        end

        function a = getDerivative(obj)
            a  = obj.alpha;
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.alpha0 = cParams.alpha0;
            obj.alpha1 = cParams.alpha1;
            obj.density = cParams.density; 
            obj.pExp = 1;        
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