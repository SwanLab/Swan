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
            s.mxMax = 0.99;
            s.myMax = 0.99;
            s.mxMin = 0.01;
            s.myMin = 0.01;
            s.xi = obj.xi;
            s.rho = obj.rho;
            sC = SuperellipseExponentBoundsComputer(s);
            obj.qMin = sC.qMin;
            obj.qMax = sC.qMax;
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
       
end