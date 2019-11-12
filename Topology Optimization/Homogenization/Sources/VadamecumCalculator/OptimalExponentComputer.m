classdef OptimalExponentComputer < handle
    
    properties (Access = public)
       qOpt 
    end
    
    properties (Access = private)
        maxStress
        mx
        my
        rho
        txi
        rhoV
        txiV
        mxMax
        myMax
        psi
        psiV
        samplePoints
        fileName
    end
    
    methods (Access = public)
        
        function obj = OptimalExponentComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            npsi = length(obj.psiV);
            npoints = length(obj.rhoV);
            obj.qOpt = zeros(npsi,npoints);
            for ipsi = 1:npsi
                obj.psi = obj.psiV(ipsi);
                disp([num2str(ipsi/npsi*100),'%'])                
                for ipoint = 1:npoints
                    obj.rho = obj.rhoV(ipoint);
                    obj.txi = obj.txiV(ipoint);
                    obj.qOpt(ipsi,ipoint) = obj.computeOptimalExponent();
                end
            end
            x.rho = obj.rhoV;
            x.txi = obj.txiV;
            x.psi = obj.psiV;
            x.q = obj.qOpt;
            save(obj.fileName,'x');
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.samplePoints = cParams.samplePoints;
            obj.computeRhoTxiAndPsiSamples();
        end
        
        function computeRhoTxiAndPsiSamples(obj)
            obj.samplePoints.compute();
            obj.rhoV = obj.samplePoints.rhoV;
            obj.txiV = obj.samplePoints.txiV;
            obj.psiV = obj.samplePoints.psiV;
        end
        
        function x = computeOptimalExponent(obj)
            f = @(q) obj.computingMaxStress(q);
            qmin = 2;
            qmax = 32;
            options = optimset('Display','iter');
            %options = optimset();            
            x = fminbnd(f,qmin,qmax,options);
        end
        
        function maxStress = computingMaxStress(obj,q)
            obj.computeMxMy(q);
            fName = 'OptimalSuperEllipse';
            outputFolder = fullfile(pwd,'Output',fName);
            gmsFile = [fullfile(outputFolder,fName),'.msh'];
            obj.createMesh(fName,outputFolder,q);
            homog = obj.createNumericalHomogenizer(fName,gmsFile);
            maxStress = obj.computeMaxStressNorm(homog);
        end
        
        function maxSnorm = computeMaxStressNorm(obj,homog)
            Ch = homog.cellVariables.Ch;
            microProblem = homog.getMicroProblem();
            stress = [cos(obj.psi) sin(obj.psi) 0]';
            strain = Ch\stress;
            microProblem.element.setVstrain(strain');
            microProblem.computeVariables;
            stresses = microProblem.variables.stress;
            sx = squeeze(stresses(:,1,:));
            sy = squeeze(stresses(:,2,:));
            sxy = squeeze(stresses(:,3,:));
            sNorm = sqrt(sx.*sx + 2*sxy.*sxy + sy.*sy);
            maxSnorm = max(sNorm);
        end
        
        function createMesh(obj,fileName,outputFolder,q)
            d = SettingsFreeFemMeshGenerator();
            d.freeFemFileName = 'SmoothRectangle';
            d.hMax  = 0.02;%0.0025;
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
            d.print   = false;
            d.iter = 0;
            nH = NumericalHomogenizerCreatorFromGmsFile(d);
            homog = nH.getHomogenizer();
        end
        
        function computeMxMy(obj,q)
            obj.computeMx(q);
            obj.computeMy(q);
        end
        
        function computeMx(obj,q)
            n = (1-obj.rho)*tan(obj.txi);
            d = obj.cFunction(q);
            obj.mx = sqrt(n/d);
        end
        
        function computeMy(obj,q)
            n = (1-obj.rho);
            d = obj.cFunction(q)*tan(obj.txi);
            obj.my = sqrt(n/d);
        end

    end
    
    methods (Access = private, Static)
        
        function c = cFunction(q)
            c = gamma(1 + 1/q)^2/gamma(1 + 2/q);
        end
        
    end
    
end



