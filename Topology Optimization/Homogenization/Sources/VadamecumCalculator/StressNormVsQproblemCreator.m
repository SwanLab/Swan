classdef StressNormVsQproblemCreator < handle
    
    properties (Access = public)
        objective
        iter
        qMin
        qMax
    end
    
    properties (Access = private)
        rho
        xi
        fileName
        phi
        hMesh
        pNorm
        print
        hasToCaptureImage
    end
    
    methods (Access = public)
        
        function obj = StressNormVsQproblemCreator(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeQbounds();
            obj.computeObjectiveFunction();
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.rho = cParams.rho;
            obj.xi  = cParams.xi;
            obj.fileName = cParams.fileName;
            obj.phi = cParams.phi;
            obj.hMesh = cParams.hMesh;
            obj.pNorm = cParams.pNorm;
        end
        
        function computeObjectiveFunction(obj)
            obj.objective = @(q) obj.computeMaximumStress(q);
        end
        
        function computeQbounds(obj)
            c = @(q) gamma(1 + 1/q)^2/gamma(1 + 2/q);
            mxMax = 0.99;
            myMax = 0.99;
            mxMin = 0.01;
            myMin = 0.01;
            obj.qMin = obj.computeQmin(c,mxMax,myMax);
            obj.qMax = obj.computeQmax(c,mxMin,myMin);
        end
        
        function qMin = computeQmin(obj,c,mxMax,myMax)
            qMin = 2;
            cqMinMx =  (1 - obj.rho)/(mxMax*tan(obj.xi));
            cqMinMy =  (1 - obj.rho)*tan(obj.xi)/myMax;
            qMinMx = fzero(@(q) c(q) - cqMinMx,16);
            qMinMy = fzero(@(q) c(q) - cqMinMy,16);
            qMin = max(qMin,max(qMinMx,qMinMy));
        end
        
        function qMax = computeQmax(obj,c,mxMin,myMin)
            qMax = 32;
            cqMaxMx =  (1 - obj.rho)/(mxMin*tan(obj.xi));
            cqMaxMy =  (1 - obj.rho)*tan(obj.xi)/(myMin);
            qMaxMx = fzero(@(q) c(q) - cqMaxMx,16);
            qMaxMy = fzero(@(q) c(q) - cqMaxMy,16);
            qMax = min(qMax,min(qMaxMx,qMaxMy));
        end
        
        function sPnorm = computeMaximumStress(obj,q)
            s.mx       = SuperEllipseParamsRelator.mx(obj.xi,obj.rho,q);
            s.my       = SuperEllipseParamsRelator.my(obj.xi,obj.rho,q);
            s.q        = q;
            s.phi      = obj.phi;
            s.pNorm    = obj.pNorm;
            s.print    = obj.print;
            s.hMesh    = obj.hMesh;
            s.iter     = obj.iter;
            s.fileName = obj.fileName;
            s.hasToCaptureImage = obj.hasToCaptureImage;
            sN = StressNormSuperEllipseComputer(s);
            sPnorm = sN.compute();
%            sN.printStress();
        end
        
    end
    
    
end