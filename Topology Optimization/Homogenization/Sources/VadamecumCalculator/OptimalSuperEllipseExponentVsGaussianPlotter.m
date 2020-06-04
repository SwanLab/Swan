classdef OptimalSuperEllipseExponentVsGaussianPlotter < handle
    
    properties (Access = private)
        samplePoints
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
        phiDV
        phiV
        iIn
        iOut
        qMin
        qMax
        qMean
    end
    
    methods (Access = public)
        
        function obj = OptimalSuperEllipseExponentVsGaussianPlotter()
            obj.init();
            for iTest = 1:4%length(obj.txiV)
                obj.rho   = obj.rhoV(iTest);
                obj.txi   = obj.txiV(iTest);
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
          %  T = obj.phiIntegrationInterval;
          %  phi = linspace((pi - obj.txi)-T,(pi - obj.txi)+T,obj.nPhi);
            phi = linspace(obj.phiMin,obj.phiMax,obj.nPhi);
        end
        
        function init(obj)
            obj.outputPath = '/home/alex/Dropbox/PaperStress/';
            %obj.phiMin = -pi;
            %obj.phiMax = pi;
            %obj.nPhi = 50;
            obj.phiMin = 0;
            obj.phiMax = pi;
            obj.nPhi = 50;%50;
            
            
            
            obj.rhoV  = [0.9,0.9,0.5,0.5];
            obj.txiV  = pi/2 - [0.1083,0.557,0.88974,1.0984];
            
            %obj.rhoV = obj.rhoV(1);
            %obj.txiV = obj.txiV(1);
            
            obj.phiIntegrationInterval = (obj.phiMax - obj.phiMin)/2;
        end
        
        function createSamplePoints(obj)
            s.type = 'FromFixedRhoAndTxi';
            s.rho0 = obj.rho;
            s.txi  = obj.txi;
            s.phi  = obj.phiDV;
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
            obj.qMax = exponentComputer.qMax;
            obj.qMin = exponentComputer.qMin;
        end
        
        function computeMeanSuperEllipseExponent(obj)
            s.phiV   = obj.phiV;
            qMeanComputer = SuperEllipseMeanExponentComputer(s);
            obj.qMean = qMeanComputer.compute(obj.qOpt,obj.txi);
        end        
       
        function plotQoptAndGaussianVsPhi(obj,itxi)
            obj.figureID = figure();
            obj.plotQvsPhi();
            obj.plotBothPvsTxiObtainPvsPhi();
            obj.addXlabel();
            obj.addTitle();
            obj.addLegend();
            %obj.print(itxi);
        end
        
        function addLegend(obj)
            l1 = '$q^*$';
            l2 = '$q_{LB}$';
            l3 = '$q_{UB}$';            
            l4 = '$q_N$';
            l5 = '$P_{\xi}$';
            l6 = '$P_{\pi - \xi}$';
            legend(l1,l2,l3,l4,l5,l6,'Interpreter','latex','Location','northeast');            
        end
        
        function plotQvsPhi(obj)
            phi(:,1) = obj.phiV;
            q(:,1)   = obj.qOpt;
            yyaxis left
            plot(phi,q,'+-','LineWidth',4);
            hold on
            plot(obj.phiV,obj.qMin,'--','LineWidth',3);    
            cs = get(gca,'colororder');            
            plot(obj.phiV,obj.qMax,'--','LineWidth',3,'Color',cs(1,:));                              
            plot([obj.phiMin,obj.phiMax],[obj.qMean,obj.qMean],'-','LineWidth',4,'Color','k');          
            axis([obj.phiMin,obj.phiMax,0,35])
        end
        
        function plotBothPvsTxiObtainPvsPhi(obj)
            xi = obj.txi;
            phi(:,1) = linspace(obj.phiMin,obj.phiMax,1000);             
            Pl = SuperEllipseMeanExponentComputer.obtainPvsPhi(phi,xi);
            Pr = SuperEllipseMeanExponentComputer.obtainPvsPhi(phi,pi - xi);            
            yyaxis right
            hold on
            plot(phi,Pl,'LineWidth',3,'LineStyle','--');
            cs = get(gca,'colororder');
            hold on
            plot(phi,Pr,'LineWidth',3,'Color',cs(1,:),'LineStyle','-.');
            ylabel('P');
        end
        
        function addXlabel(obj)
            xlabel('$\phi$','Interpreter','latex');
             set(gca,'xtick',[obj.phiMin:pi/8:obj.phiMax])
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2',...
                                    '5\pi/8','3\pi/4','7\pi/8','\pi'})            
        end
        
        function addTitle(obj)
            rhoT = ['\rho = ',num2str(obj.rho)];
            txiT = ['\xi = ',num2str(obj.txi)];
            tit =  ['$',rhoT,'\quad',txiT,'$'];
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