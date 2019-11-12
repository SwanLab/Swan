classdef OptimalSuperEllipseExponentVsGaussianPlotter < handle
    
    properties (Access = private)
        samplePoints
        gaussian
        qOpt
        psi
        psiV
        icase
        rho0
        figureID
        outputPath        
    end
    
    methods (Access = public)
        
        function obj = OptimalSuperEllipseExponentVsGaussianPlotter()
            obj.init();
            for iphi = 1:length(obj.psiV)
                obj.icase = iphi;
                obj.psi = obj.psiV(iphi);
                obj.createSamplePoints();
                obj.computeOptimalSuperEllipseExponent();
                obj.plotQoptAndGaussianVsTxi(); 
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.outputPath = '/home/alex/Dropbox/PaperStress/';            
            obj.gaussian = gaussianFunction();
            obj.psiV = linspace(0,pi/2,6);
            obj.rho0 = 0.8;
        end
        
        function createSamplePoints(obj)
            s.type = 'FromFixedRho';
            s.rho0 = obj.rho0;
            s.psi  = obj.psi;
            sample = SamplePointsCreatorForOptimalExponentComputer.create(s);
            obj.samplePoints = sample;
        end
        
        function computeOptimalSuperEllipseExponent(obj)
            s.samplePoints = obj.samplePoints;
            s.fileName = 'OptimalSuperEllipseExponentDataFromFixedRho';
            exponentComputer = OptimalExponentComputer(s);
            exponentComputer.compute();
            obj.qOpt = exponentComputer.qOpt;
        end
        
        function plotQoptAndGaussianVsTxi(obj)
            obj.figureID = figure();
            obj.plotQvsTxi();
            obj.plotPvsTxiobtainPvsTxi();
            obj.addXlabel();
            obj.print()
        end               
        
        function plotQvsTxi(obj)
            obj.samplePoints.compute();
            txi(:,1) = obj.samplePoints.txiV;
            q(:,1)   = obj.qOpt;   
            yyaxis left            
            plot(txi,q,'+-','LineWidth',4);
            ylabel('q');            
        end
        
        function plotPvsTxiobtainPvsTxi(obj)
            npoints = 100;
            phiValue = obj.psi;
            txiMin = obj.samplePoints.txiMin;
            txiMax = obj.samplePoints.txiMax;            
            txi =  linspace(txiMin,txiMax,npoints);
            P = obj.gaussian(txi,phiValue);  
            yyaxis right
            plot(txi,P,'LineWidth',4)
            ylabel('P');            
        end
        
        function addXlabel(obj)
            xlabel('$\xi$','Interpreter','latex');            
            set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})              
        end
        
        function print(obj)
            fp = contourPrinter(obj.figureID);            
            filePath = [obj.outputPath,'QoptAndGaussianVsTxi',num2str(obj.icase)];
            fp.print(filePath);             
        end
        
    end
    
end