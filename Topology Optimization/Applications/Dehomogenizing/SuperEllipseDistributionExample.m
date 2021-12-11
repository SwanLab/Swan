classdef SuperEllipseDistributionExample < handle
    
    properties (Access = public)
      m1
      m2
      q
    end
    
    properties (Access = private)
      
    end
    
    properties (Access = private)
       coord 
       mMin
       mMax
       qMin
       qMax
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseDistributionExample(cParams)
            obj.init(cParams)
        end
        
        function computeParameters(obj)
            obj.computeM1M2();
            obj.computeSmoothExponent();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.coord    = cParams.coord;
           obj.mMin     = cParams.mMin;
           obj.mMax     = cParams.mMax;
           obj.qMin     = cParams.qMin;
           obj.qMax     = cParams.qMax;
        end
        
        function computeM1M2(obj)
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);
            p = 16;
            r = (x1.^p + x2.^p).^(1/p);
            r = r.^(1/p);
            obj.m1 = obj.createLinearFunction(r,obj.mMin,obj.mMax);
            obj.m2 = obj.createLinearFunction(r,obj.mMin,obj.mMax);
        end
        
        function computeSmoothExponent(obj)
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);
            xM = max(x1,x2);
            obj.q = obj.createLinearFunction(xM,obj.qMin,obj.qMax);
        end
        
    end
    
    methods (Access = private, Static)
       
        function f = createLinearFunction(x,fMin,fMax)
            xmin = min(x);
            xmax = max(x);
            m = fMin +(fMax-fMin)*(x-xmin)/(xmax-xmin);
            f = (m);
        end
        
    end
    
end