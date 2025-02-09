classdef OptimalSuperEllipseExponentVsGaussianPlotter < handle
    
    properties (Access = private)
       figureID      
       phiMin
       phiMax        
    end    
    
    properties (Access = private)
       phiV
       qOpt
       qMin
       qMax
       qMean
       xi
       rho
       outputPath
    end
    
    methods (Access = public)
        
        function obj = OptimalSuperEllipseExponentVsGaussianPlotter(cParams)
            obj.init(cParams)            
        end                
        
         function plot(obj,iPlot)
            obj.figureID = figure();
            obj.plotQvsPhi();
            obj.plotBothPvsTxiObtainPvsPhi();
            obj.addXlabel();
            obj.addTitle();
            obj.addLegend();
            obj.print(iPlot);
        end        
        
    end
    
    methods (Access = private)
       
        function init(obj,s)
            obj.phiV   = s.phiV;
            obj.qOpt   = s.qOpt;
            obj.qMin   = s.qMin;
            obj.qMax   = s.qMax;
            obj.qMean  = s.qMean;
            obj.xi     = s.xi;
            obj.rho    = s.rho; 
            obj.phiMin = min(obj.phiV);            
            obj.phiMax = max(obj.phiV);
            obj.outputPath = '/home/alex/Dropbox/PaperStress/OptimalSmoothingPonderated/';
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
            plot([obj.phiMin,obj.phiMax],[obj.qMin,obj.qMin],'--','LineWidth',3);    
            cs = get(gca,'colororder');            
            plot([obj.phiMin,obj.phiMax],[obj.qMax,obj.qMax],'--','LineWidth',3,'Color',cs(1,:));                              
            plot([obj.phiMin,obj.phiMax],[obj.qMean,obj.qMean],'-','LineWidth',4,'Color','k');          
            axis([obj.phiMin,obj.phiMax,0,35])
        end
        
        function plotBothPvsTxiObtainPvsPhi(obj)
            phi(:,1) = linspace(obj.phiMin,obj.phiMax,1000);             
            Pl = SuperEllipseMeanExponentComputer.obtainPvsPhi(phi,obj.xi);
            Pr = SuperEllipseMeanExponentComputer.obtainPvsPhi(phi,pi - obj.xi);            
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
            xValues = obj.phiMin:pi/8:obj.phiMax;
            set(gca,'xtick',xValues)
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2',...
                                    '5\pi/8','3\pi/4','7\pi/8','\pi'})            
        end
        
        function addTitle(obj)
            rhoT = ['\rho = ',num2str(obj.rho)];
            txiT = ['\xi = ',num2str(obj.xi*180/pi)];
            tit =  ['$',rhoT,'\quad',txiT,'$'];
            title(tit,'Interpreter','latex');
        end
        
        function print(obj,iPlot)
            fp = contourPrinter(obj.figureID);
            filePath = [obj.outputPath,'QoptAndGaussianVsPhi',num2str(iPlot)];
            fp.print(filePath);
        end        
        
        
    end
    
end