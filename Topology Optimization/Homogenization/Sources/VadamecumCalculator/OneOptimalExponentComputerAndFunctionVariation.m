classdef OneOptimalExponentComputerAndFunctionVariation < handle
    
    properties (Access = public)
        qOptIter
        fOptIter
        fValues
        qValues
    end
    
    properties (Access = private)
        qMin
        qMax
    end
    
    properties (Access = private)
        rho
        txi
        phi
        pNorm
        print
        hasToCaptureImage
        fileName
        hMesh
    end
    
    methods (Access = public)
        
        function obj = OneOptimalExponentComputerAndFunctionVariation(cParams)
            obj.init(cParams);
            obj.computeQbounds();
        end
        
        function compute(obj)
            obj.computeOptimalExponent();
            obj.printOptimalMicroStructure();
            obj.computeStressNormRelationWithQ();
        end
        
        function printOptimalMicroStructure(obj)
            obj.print = true;
            obj.hasToCaptureImage = true;
            qOpt = obj.qOptIter(end);
            obj.computingMaxStress(qOpt);
        end
        
        function computeOptimalExponent(obj)
            f = @(q) obj.computingMaxStress(q);
            xIter = [];
            fIter = [];
            options = optimset('Display','iter','OutputFcn', @myoutput);
            options.MaxIter = 200;
            options.TolX = 1e-5;
            x = fminbnd(f,obj.qMin,obj.qMax,options);
            function stop = myoutput(x,optimvalues,state)
                stop = false;
                if isequal(state,'iter')
                    xIter = [xIter; x];
                    fIter = [fIter;optimvalues.fval];
                end
            end
            obj.qOptIter = xIter;
            obj.fOptIter = fIter;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.rho   = cParams.rho;
            obj.txi   = cParams.txi;
            obj.phi   = cParams.phi;
            obj.pNorm = cParams.pNorm;
            obj.hMesh = cParams.hMesh;
            obj.print = false;
            obj.hasToCaptureImage = false;
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
        
        function computeStressNormRelationWithQ(obj)
            nValues = 50;
            qVal = linspace(obj.qMin,obj.qMax,nValues);
            fVal = zeros(nValues,1);
            for ivalues = 2:nValues
                qV = qVal(ivalues);
                fVal(ivalues) = obj.computingMaxStress(qV);
            end
            obj.qValues = qVal;
            obj.fValues = fVal;
        end
        
        function qMin = computeQmin(obj,c,mxMax,myMax)
            qMin = 2;
            cqMinMx =  (1 - obj.rho)/(mxMax*tan(obj.txi));
            cqMinMy =  (1 - obj.rho)*tan(obj.txi)/myMax;
            qMinMx = fzero(@(q) c(q) - cqMinMx,16);
            qMinMy = fzero(@(q) c(q) - cqMinMy,16);
            qMin = max(qMin,max(qMinMx,qMinMy));
        end
        
        function qMax = computeQmax(obj,c,mxMin,myMin)
            qMax = 32;
            cqMaxMx =  (1 - obj.rho)/(mxMin*tan(obj.txi));
            cqMaxMy =  (1 - obj.rho)*tan(obj.txi)/(myMin);
            qMaxMx = fzero(@(q) c(q) - cqMaxMx,16);
            qMaxMy = fzero(@(q) c(q) - cqMaxMy,16);
            qMax = min(qMax,min(qMaxMx,qMaxMy));
        end
        
        function sPnorm = computingMaxStress(obj,q)
            s.mx       = obj.computeMx(q);
            s.my       = obj.computeMy(q);
            s.q        = q;
            s.phi      = obj.phi;
            s.pNorm    = obj.pNorm;
            s.print    = obj.print;
            s.hMesh    = obj.hMesh;
            s.fileName = [obj.fileName,'Q',strrep(num2str(q),'.','_')];
            s.hasToCaptureImage = obj.hasToCaptureImage;
            sN = StressNormSuperEllipseComputer(s);
            sPnorm = sN.compute();
        end
        
        function mx = computeMx(obj,q)
            mx = SuperEllipseParamsRelator.mx(obj.txi,obj.rho,q);
        end
        
        function my = computeMy(obj,q)
            my = SuperEllipseParamsRelator.my(obj.txi,obj.rho,q);
        end
        
    end
    
end


