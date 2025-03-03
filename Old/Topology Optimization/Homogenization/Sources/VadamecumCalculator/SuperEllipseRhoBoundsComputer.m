classdef SuperEllipseRhoBoundsComputer < handle
    
    properties (Access = private)
       q
       txi
       mxMax
       myMax
       mxMin
       myMin
       superEllipse
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
            obj.superEllipse = SuperEllipseParamsRelator;
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
            rho = obj.superEllipse.rhoFromMxAndTxi(mx,obj.txi,obj.q);
        end

        function rho = rhoFromMy(obj,my)
            rho = obj.superEllipse.rhoFromMyAndTxi(my,obj.txi,obj.q);
        end                
        
    end
    
    
end