classdef SuperEllipseRhoBoundsComputer < handle
    
    properties (Access = private)
       q
       txi
       mxMax
       myMax
       mxMin
       myMin
    end
    
    methods (Access = public)
       
        function obj = SuperEllipseRhoBoundsComputer(cParams)
            obj.init(cParams);
        end
        
        function [rhoMin,rhoMax] = compute(obj)
            rhoMin = obj.computeRhoMin();
            rhoMax = obj.computeRhoMax();            
        end
            
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.q = cParams.q;
            obj.txi = cParams.txi;
            obj.mxMin = cParams.mxMin;
            obj.myMin = cParams.myMin;  
            obj.mxMax = cParams.mxMax;            
            obj.myMax = cParams.myMax;                        
        end

        function rhoMin = computeRhoMin(obj)         
            rhoMinMx = obj.rhoFromMx(obj.mxMax);
            rhoMinMy = obj.rhoFromMy(obj.myMax);
            rhoMin = max(rhoMinMx,rhoMinMy);
        end
        
        function rhoMax = computeRhoMax(obj)         
            rhoMaxMx = obj.rhoFromMx(obj.mxMin);
            rhoMaxMy = obj.rhoFromMy(obj.myMin);
            rhoMax = min(rhoMaxMx,rhoMaxMy);
        end
        
        function rho = rhoFromMx(obj,mx)
          rho = 1 - obj.c(obj.q)*mx^2*tan(obj.txi);
        end

        function rho = rhoFromMy(obj,my)
          rho = 1 - obj.c(obj.q)*my^2./tan(obj.txi);
        end                
        
    end
    
    methods (Access = private, Static)
        
        function c = c(q)
            c = gamma(1 + 1/q)^2/gamma(1 + 2/q);            
        end        
        
    end
    
    
end