classdef SomeOptimalSuperEllipseExponentVsStress < handle
    
    properties (Access = private)
        samplePoints
        qOpt
        phiMin
        phiMax
        nPhi
        rhoV
        xiV
        rho
        xi
        outputPath
        phiDV
        phiV
        iIn
        iOut
        qMin
        qMax
        qMean
    end
    
    methods (Access = public)
        
        function obj = SomeOptimalSuperEllipseExponentVsStress()
            obj.init();
            for iTest = 1:length(obj.xiV)
                obj.rho   = obj.rhoV(iTest);
                obj.xi   = obj.xiV(iTest);
                obj.phiV  = obj.createPhi();
                obj.phiDV = obj.createPhiInInterval();       
                obj.createSamplePoints();
                obj.computeOptimalSuperEllipseExponent();
                obj.computeMeanSuperEllipseExponent();            
                obj.plotQoptAndGaussianVsPhi(iTest);                
            end
        end
        
    end    
    
    methods (Access = private)
        
        function init(obj)
            obj.outputPath = '/home/alex/Dropbox/PaperStress/';
            obj.phiMin = 0;
            obj.phiMax = pi;
            obj.nPhi = 3;%50;               
            obj.rhoV  = [0.9,0.9,0.5,0.5];
            obj.xiV  = pi/2 - [0.1083,0.557,0.88974,1.0984];            
        end        
        
        function itIs = isInInterval(obj,phi)
            itIs = phi>= obj.phiMin & phi <= obj.phiMax;
        end
        
        function phiD = createPhiInInterval(obj)
            phi = obj.phiV;            
            iI  = obj.isInInterval(phi);
            iO  = ~iI;
            phiD(iI) = phi(iI);
            phiD(iO) = sign(phi(iO)).*(phi(iO)-pi/2);            
        end
        
        function phi = createPhi(obj)
            phi = linspace(obj.phiMin,obj.phiMax,obj.nPhi);
        end
         
        function createSamplePoints(obj)
            s.type = 'FromFixedRhoAndTxi';
            s.rho0 = obj.rho;
            s.txi  = obj.xi;
            s.phi  = obj.phiV;
            sample = SamplePointsCreatorForOptimalExponentComputer.create(s);
            sample.compute();
            obj.samplePoints = sample;
        end
        
        function computeOptimalSuperEllipseExponent(obj)
            s.samplePoints = obj.samplePoints;
            rhoT = strrep(num2str(obj.rho),'.','_');
            txiT = strrep(num2str(round(obj.xi,3)),'.','_');
            fN = ['AveragingSuperEllipseRho',rhoT,'Txi',txiT];
            s.fileName = fN;
            exponentComputer = OptimalExponentComputer(s);
            exponentComputer.compute();
            obj.qOpt = exponentComputer.qOpt;
            obj.qMax = exponentComputer.qMax;
            obj.qMin = exponentComputer.qMin;
        end
        
        function computeMeanSuperEllipseExponent(obj)
            s.phiV = obj.phiV;
            qMeanC    = SuperEllipseMeanAndDesvExponentComputer(s);
            obj.qMean = qMeanC.computeMean(obj.qOpt,obj.xi);
        end        
        
        function plotQoptAndGaussianVsPhi(obj,itxi)
            s.phiV = obj.phiV;
            s.phiMin = obj.phiMin;
            s.phiMax = obj.phiMax;
            s.qOpt = obj.qOpt;
            s.qMin = obj.qMin;
            s.qMax = obj.qMax;
            s.qMean = obj.qMean;
            s.xi = obj.xi;
            s.rho = obj.rho;
            p = OptimalSuperEllipseExponentVsGaussianPlotter(s);
            p.plot(itxi);
        end
        
    end
    
end