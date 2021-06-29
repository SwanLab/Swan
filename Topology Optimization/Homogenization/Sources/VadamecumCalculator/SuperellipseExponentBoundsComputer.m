classdef SuperellipseExponentBoundsComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        qMin
        qMax
    end
    
    properties (Access = private)
        mxMax
        myMax
        mxMin
        myMin
        xi
        rho
        qMin0
        qMax0        
    end
    
    methods (Access = public)
        
        function obj = SuperellipseExponentBoundsComputer(cParams)
            obj.init(cParams)
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mxMax = cParams.mxMax;
            obj.myMax = cParams.myMax;
            obj.mxMin = cParams.mxMin;
            obj.myMin = cParams.myMin;
            obj.rho   = cParams.rho;
            obj.xi    = cParams.xi;
            obj.qMin0 = 2;
            obj.qMax0 = 32;
        end
        
        function compute(obj)
            c = @(q) gamma(1 + 1/q)^2/gamma(1 + 2/q);
            obj.qMin = obj.computeQmin(c,obj.mxMax,obj.myMax);
            obj.qMax = obj.computeQmax(c,obj.mxMin,obj.myMin);
        end
        
        function qMin = computeQmin(obj,c,mxMax,myMax)
            cqMinMx =  (1 - obj.rho)/(mxMax*tan(obj.xi));
            cqMinMy =  (1 - obj.rho)*tan(obj.xi)/myMax;
            qMinMx = obj.obtainQsolution(c,cqMinMx);
            qMinMy = obj.obtainQsolution(c,cqMinMy);
            qMin = max(obj.qMin0,max(qMinMx,qMinMy));
        end
        
        function qMax = computeQmax(obj,c,mxMin,myMin)
            cqMaxMx =  (1 - obj.rho)/(mxMin*tan(obj.xi));
            cqMaxMy =  (1 - obj.rho)*tan(obj.xi)/(myMin);
            qMaxMx = obj.obtainQsolution(c,cqMaxMx);
            qMaxMy = obj.obtainQsolution(c,cqMaxMy);
            qMax = min(obj.qMax0,min(qMaxMx,qMaxMy));
        end
        
    end
    
   methods (Access = private, Static)
       
        function q = obtainQsolution(c,cV)
            problem.x0 = 16;
            problem.options = optimset('Display','off');    
            problem.objective = @(q) c(q) - cV;
            problem.solver = 'fzero';
            q = fzero(problem);                    
        end
        
    end    
       
end