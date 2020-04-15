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
        mx
        my
    end
    
    properties (Access = private)
        rho
        txi
        psi        
        pNorm
        print
        hasToCaptureImage    
        fileName
        hMesh
    end
    
    methods (Access = public)
        
        function obj = OneOptimalExponentComputerAndFunctionVariation(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
           obj.computeQbounds();
           obj.computeOptimalExponent();
           obj.computeStressNormRelationWithQ();
           obj.printOptimalMicroStructure();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.rho   = cParams.rho;
            obj.txi   = cParams.txi;
            obj.psi   = cParams.psi;   
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
        
        function computeOptimalExponent(obj)
            f = @(q) obj.computingMaxStress(q);
            xIter = [];
            fIter = [];
            options = optimset('Display','iter','OutputFcn', @myoutput);
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
        
        function maxStress = computingMaxStress(obj,q)
            obj.computeMxMy(q);
            fName = ['OptimalSuperEllipse',obj.fileName];
            outputFolder = fullfile(pwd,'Output',fName);
            gmsFile = [fullfile(outputFolder,fName),'.msh'];
            obj.createMesh(fName,outputFolder,q);
            homog = obj.createNumericalHomogenizer(fName,gmsFile);
            maxStress = obj.computeMaxStressNorm(homog);
        end
        
        function sPnorm = computeMaxStressNorm(obj,homog)
            Ch = homog.cellVariables.Ch;
            microProblem = homog.getMicroProblem();
            stress = [cos(obj.psi) sin(obj.psi) 0]';
            strain = Ch\stress;
            microProblem.element.setVstrain(strain');
            microProblem.computeVariables;
            stresses = microProblem.variables.stress;
            m = microProblem.mesh;
            q = Quadrature.set(m.geometryType);
            q.computeQuadrature('CONSTANT');
            dV = m.computeDvolume(q);
            sx  = squeeze(stresses(:,1,:));
            sy  = squeeze(stresses(:,2,:));
            sxy = squeeze(stresses(:,3,:));
            sNorm2 = sqrt(sx.*sx + 2*sxy.*sxy + sy.*sy);
            if isequal(obj.pNorm,'max')            
                sPnorm = max(sNorm2);                
            else
                p = obj.pNorm;
                int = sNorm2.^p;
                sPnorm = sum(int(:).*dV(:))^(1/p);
            end
        end
        
        function createMesh(obj,fileName,outputFolder,q)
            d = SettingsFreeFemMeshGenerator();
            d.freeFemFileName = 'SmoothRectangle';
            d.hMax  = obj.hMesh;%0.002;%0.0025;
            d.mxV             = obj.mx;
            d.myV             = obj.my;
            d.fileName        = fileName;
            d.printingDir     = outputFolder;
            d.qNorm           = q;
            fG = FreeFemMeshGenerator(d);
            fG.generate();
        end
        
        function homog = createNumericalHomogenizer(obj,fileName,gmsFile)
            d.gmsFile = gmsFile;
            d.outFile = fileName;
            d.print   = obj.print;
            d.iter = 0;
            d.hasToCaptureImage = obj.hasToCaptureImage;
            nH = NumericalHomogenizerCreatorFromGmsFile(d);
            homog = nH.getHomogenizer();
        end
        
        function printOptimalMicroStructure(obj)
            obj.print = true;
            obj.hasToCaptureImage = true;
            qOpt = obj.qOptIter(end);
            obj.computingMaxStress(qOpt);
        end            
        
        function computeMxMy(obj,q)
            obj.computeMx(q);
            obj.computeMy(q);
        end
        
        function computeMx(obj,q)
            obj.mx = SuperEllipseParamsRelator.mx(obj.txi,obj.rho,q);            
        end
        
        function computeMy(obj,q)
            obj.my = SuperEllipseParamsRelator.my(obj.txi,obj.rho,q);
        end

    end    
    
end



