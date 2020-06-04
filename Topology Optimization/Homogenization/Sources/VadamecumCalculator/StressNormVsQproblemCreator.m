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
        stressNormParameters
    end
    
    methods (Access = public)
        
        function obj = StressNormVsQproblemCreator(cParams)
            obj.init(cParams);
            obj.computeStressNormParameters();
        end
        
        function compute(obj)
            obj.computeQbounds();
            obj.computeObjectiveFunction();
        end
        
        function var = computeCellVariables(obj,q)
            s = obj.stressNormParameters;
            s.mx = SuperEllipseParamsRelator.mx(obj.xi,obj.rho,q);
            s.my = SuperEllipseParamsRelator.my(obj.xi,obj.rho,q);
            s.q  = q;
            sN = StressNormSuperEllipseComputer(s);       
            var = sN.computeCellVariables();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.rho      = cParams.rho;
            obj.xi       = cParams.xi;
            obj.fileName = cParams.fileName;
            obj.phi      = cParams.phi;
            obj.hMesh    = cParams.hMesh;
            obj.pNorm    = cParams.pNorm;
            obj.print    = cParams.print;
            obj.hasToCaptureImage = cParams.hasToCaptureImage;           
        end

        function computeStressNormParameters(obj)
            s.phi      = obj.phi;
            s.pNorm    = obj.pNorm;
            s.print    = obj.print;
            s.hMesh    = obj.hMesh;
            s.iter     = obj.iter;
            s.fileName = obj.fileName;            
            s.hasToCaptureImage = obj.hasToCaptureImage;            
            obj.stressNormParameters = s;
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
            qMinMx = obj.obtainQsolution(c,cqMinMx);
            qMinMy = obj.obtainQsolution(c,cqMinMy);
            qMin = max(qMin,max(qMinMx,qMinMy));
        end
        
        function qMax = computeQmax(obj,c,mxMin,myMin)
            qMax = 32;
            cqMaxMx =  (1 - obj.rho)/(mxMin*tan(obj.xi));
            cqMaxMy =  (1 - obj.rho)*tan(obj.xi)/(myMin);
            qMaxMx = obj.obtainQsolution(c,cqMaxMx);
            qMaxMy = obj.obtainQsolution(c,cqMaxMy);
            qMax = min(qMax,min(qMaxMx,qMaxMy));
        end
        
        function sPnorm = computeMaximumStress(obj,q)
            s = obj.stressNormParameters;
            s.mx = SuperEllipseParamsRelator.mx(obj.xi,obj.rho,q);
            s.my = SuperEllipseParamsRelator.my(obj.xi,obj.rho,q);
            s.q  = q;
            sN = StressNormSuperEllipseComputer(s);
            sPnorm = sN.compute();
            sN.printImage();
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