classdef SuperEllipseDistributionExample < handle
    
    properties (Access = public)
      theta
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
       thetaMin
       thetaMax
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseDistributionExample(cParams)
            obj.init(cParams)            
        end
        
        function computeParameters(obj)
            obj.computeM1M2();
            obj.computeTheta();
            obj.computeSmoothExponent();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.coord = cParams.coord;
           obj.mMin = cParams.mMin;
           obj.mMax = cParams.mMax;
           obj.qMin = cParams.qMin;
           obj.qMax = cParams.qMax;
           obj.thetaMin = cParams.thetaMin;
           obj.thetaMax = cParams.thetaMax;
        end
        
        function computeM1M2(obj)
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);
            obj.m1 = obj.createLinearFunction(x1,obj.mMin,obj.mMax);
            obj.m2 = obj.createLinearFunction(x2,obj.mMin,obj.mMax);
        end
        
        function computeTheta(obj)                      
            x1 = obj.coord(:,1);
            x2 = obj.coord(:,2);            
            xM = x1;%max(x1,x2);            
            obj.theta = obj.createLinearFunction(xM,obj.thetaMin,obj.thetaMax);
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
            m = fMax +(fMin-fMax)*(x-xmin)/(xmax-xmin);
            f = (m);
        end
        
    end    
    
end