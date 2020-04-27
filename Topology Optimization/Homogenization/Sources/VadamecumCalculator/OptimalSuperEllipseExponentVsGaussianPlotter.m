classdef OptimalSuperEllipseExponentVsGaussianPlotter < handle
    
    properties (Access = private)
        samplePoints
        gaussian
        qOpt
        phiMin
        phiMax
        nPhi
        rhoV
        txiV
        rho
        txi
        figureID
        outputPath
        phiIntegrationInterval
    end
    
    methods (Access = public)
        
        function obj = OptimalSuperEllipseExponentVsGaussianPlotter()
            obj.init();
            for itxi = 1:length(obj.txiV)
                obj.rho = obj.rhoV(itxi);
                obj.txi = obj.txiV(itxi);
                obj.createSamplePoints();
                obj.computeOptimalSuperEllipseExponent();
                obj.plotQoptAndGaussianVsPhi(itxi);
            end
        end
        
    end
    methods (Access = private)
        
        function init(obj)
            obj.outputPath = '/home/alex/Dropbox/PaperStress/';
            obj.gaussian = gaussianFunction();
            obj.phiMin = -pi;%0;
            obj.phiMax = pi;%pi/2;
            obj.nPhi = 50;
            %
            %             obj.phiMin = 0-pi/6;%0;
            %             obj.phiMax = 0+pi/6;%pi/2;
            %             obj.nPhi   = 2;
            %
            obj.rhoV  = [0.9,0.9,0.5,0.5];
            obj.txiV  = pi/2 - [0.1083,0.557,0.88974,1.0984];
            
            obj.rhoV = obj.rhoV(1);
            obj.txiV = obj.txiV(1);
            
            obj.phiIntegrationInterval = pi/4;
        end
        
        function createSamplePoints(obj)
            s.type = 'FromFixedRhoAndTxi';
            s.rho0 = obj.rho;
            s.txi  = obj.txi;
            s.phi  = linspace(obj.phiMin,obj.phiMax,obj.nPhi);
            sample = SamplePointsCreatorForOptimalExponentComputer.create(s);
            sample.compute();
            obj.samplePoints = sample;
        end
        
        function computeOptimalSuperEllipseExponent(obj)
            s.samplePoints = obj.samplePoints;
            rhoT = strrep(num2str(obj.rho),'.','_');
            txiT = strrep(num2str(round(obj.txi,3)),'.','_');
            fN = ['AveragingSuperEllipseRho',rhoT,'Txi',txiT];
            s.fileName = fN;
            exponentComputer = OptimalExponentComputer(s);
            exponentComputer.compute();
            obj.qOpt = exponentComputer.qOpt;
        end
        
        function plotQoptAndGaussianVsPhi(obj,itxi)
            obj.figureID = figure();
            obj.plotQvsPhi();
            obj.plotPvsTxiObtainPvsPhi();
            obj.addXlabel();
            obj.addTitle();
            %  obj.print(itxi);
        end
        
        function plotQvsPhi(obj)
            sE = SuperEllipseParamsRelator;
            qMax = 32;
            mxMax = 0.99;
            myMax = 0.99;
            txiMax = sE.txiFromMxAndRho(mxMax,obj.rho,qMax);
            txiMin = sE.txiFromMyAndRho(myMax,obj.rho,qMax);
            phi(:,1) = obj.samplePoints.phiV;
            q(:,1)   = obj.qOpt();
            yyaxis left
            plot(phi,q,'+-','LineWidth',4);
            ylabel('q');
            %set(gca,'xtick',[0:pi/8:pi/2])
            %set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})
            set(gca,'xtick',[obj.phiMin:pi/8:obj.phiMax])
            set(gca,'xticklabels',{'\pi','-7\pi/8','-3\pi/4','-5\pi/8',...
                '-\pi/2','-3\pi/8','-\pi/4','-\pi/8',...
                '0','\pi/8','\pi/4','3\pi/8','\pi/2',...
                '5\pi/8','3\pi/4','7\pi/8','\pi'})
            
            T = obj.phiIntegrationInterval/2;
            pMin = min(txiMin - T/2,obj.phiMin);
            pMax = max(txiMax + T/2,obj.phiMax);
            axis([pMin,pMax,2,32])
        end
        
        function plotPvsTxiObtainPvsPhi(obj)
            T   = obj.phiIntegrationInterval;
            npoints = 100;
            txiG = pi/2 - obj.txi;
            phiV = linspace(txiG-T,txiG+T,npoints);
            P = obj.gaussian(txiG,phiV);
            yyaxis right
            plot(phiV,P,'LineWidth',4)
            ylabel('P');
        end
        
        function addXlabel(obj)
            xlabel('$\phi$','Interpreter','latex');
            set(gca,'xtick',[0:pi/8:pi/2])
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})
        end
        
        function addTitle(obj)
            rhoT = ['\rho = ',num2str(obj.rho)];
            txiT = ['\xi = ',num2str(obj.txi)];
            tit =  ['$',rhoT,'\quad',txiT,'$'];%,', \quad','\xi = ',num2str(obj.txi),'$'];
            title(tit,'Interpreter','latex');
        end
        
        function print(obj,itxi)
            fp = contourPrinter(obj.figureID);
            rhoStr = strrep(num2str(obj.rho),'.','_');
            filePath = [obj.outputPath,'QoptAndGaussianVsTxi',num2str(itxi),rhoStr];
            fp.print(filePath);
        end
        
    end
    
end