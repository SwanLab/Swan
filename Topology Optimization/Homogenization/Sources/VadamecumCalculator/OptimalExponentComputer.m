classdef OptimalExponentComputer < handle
    
    properties (Access = private)
        maxStress
        mxV
        myV
        mx
        my
        rho
        txi
        rhoV
        txiV
        mxMax
        myMax
        qOpt
        psi
        psiV
    end
    
    methods (Access = public)
        
        function obj = OptimalExponentComputer()
            obj.init();
            obj.compute();
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
            save('OptimalSuperEllipseExponentData','x');
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.mxMax = 0.99;
            obj.myMax = 0.99;
            obj.mxV = linspace(0.01,obj.mxMax,10);
            obj.myV = linspace(0.01,obj.myMax,10);
            obj.psiV = linspace(pi/4,pi/4,10);
            obj.computeRhoAndTxiSamples();
        end
        
        function computeRhoAndTxiSamples(obj)
            nmx = length(obj.mxV);
            nmy = length(obj.myV);
            q = 10^6;
            for imx = 1:nmx
                for imy = 1:nmy
                    mx = obj.mxV(imx);
                    my = obj.myV(imy);
                    c = obj.cFunction(q);
                    index = nmx*(imy - 1) + imx;
                    obj.rhoV(index) = 1 - c*mx*my;
                    obj.txiV(index) = atan(mx/my);
                end
            end
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
            fileName = 'OptimalSuperEllipse';
            outputFolder = fullfile(pwd,'Output',fileName);
            gmsFile = [fullfile(outputFolder,fileName),'.msh'];
            obj.createMesh(fileName,outputFolder,q);
            homog = obj.createNumericalHomogenizer(fileName,gmsFile);
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



