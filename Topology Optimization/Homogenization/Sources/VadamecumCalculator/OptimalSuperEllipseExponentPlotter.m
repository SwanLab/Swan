classdef OptimalSuperEllipseExponentPlotter < handle
    
    properties (Access = private)
        outputPath
        optimalSuperEllipse
    end
    
    methods (Access = public)
        
        function obj = OptimalSuperEllipseExponentPlotter()
            obj.outputPath = '/home/alex/Dropbox/PaperStress/';
            obj.computePonderatedOptimalSuperEllipse();
            obj.saveQmeanMxMy();
            obj.plotQmeanTxiRho();
            obj.plotQmeanMxMy();  
        end
        
    end
    
    methods (Access = private)
        
        function computePonderatedOptimalSuperEllipse(obj)
            obj.optimalSuperEllipse = PonderatedOptimalSuperEllipseComputer();
            obj.optimalSuperEllipse.compute();
        end
                
        function plotQmeanTxiRho(obj)  
            obj.plotQmeanTxiRhoContour();
            obj.plotQmeanTxiRhoTrisurf();                        
        end        
        
        function plotQmeanTxiRhoContour(obj)
            x = obj.optimalSuperEllipse.txiV;
            y = obj.optimalSuperEllipse.rhoV;
            z = obj.optimalSuperEllipse.qMean;            
            tri = delaunay(x,y);
            f = figure();
            ncolors = 50;
            tricontour(tri,x,y,z,ncolors);            
            colorbar
            hold on
            plot(x,y,'+');
            ylim([0 1])
            xlabel('$\xi$','Interpreter','latex');
            ylabel('\rho');
            tN = '\textrm{Optimal SuperEllipse exponent}';
            title(['$',tN,'$'],'interpreter','latex')
            set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})  
            fp = contourPrinter(f);            
            filePath = [obj.outputPath,'QmeanTxiRhoContour'];
            fp.print(filePath);            
        end
        
        function saveQmeanMxMy(obj)
           n = 20;
           xmin = 0.011;0.009;
           xmax = 0.989;0.991;
           ymin = 0.011;0.009;
           ymax = 0.989;0.991;            
           [mx,my,q] = obj.interpolateQmeanMxMy(n,xmin,xmax,ymin,ymax);           
           d.mx = mx;
           d.my = my;
           d.q = q;
           fN = 'OptimalSuperEllipseExponent';
           pD = 'Topology Optimization/Vademecums';
           file2SaveName = [pD,'/',fN,'.mat'];
           save(file2SaveName,'d');            
        end
        
        function plotQmeanTxiRhoTrisurf(obj)
            x = obj.optimalSuperEllipse.txiV;
            y = obj.optimalSuperEllipse.rhoV;
            z = obj.optimalSuperEllipse.qMean; 
            tri = delaunay(x,y);            
            f = figure();
            trisurf(tri,x,y,z);
            hold on
            plot3(x,y,z,'+')
            view(-161,30)            
            xlabel('$\xi$','Interpreter','latex');
            ylabel('\rho');
            zlabel('q');
            tN = '\textrm{Optimal SuperEllipse exponent}';
            title(['$',tN,'$'],'interpreter','latex')
            set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})  
            fp = surfPrinter(f);            
            filePath = [obj.outputPath,'QmeanTxiRhoSurf'];
            fp.print(filePath);            
        end
        
        function [xq,yq,vq] = interpolateQmeanMxMy(obj,n,xmin,xmax,ymin,ymax)
            x = obj.optimalSuperEllipse.mxV;
            y = obj.optimalSuperEllipse.myV;
            v = obj.optimalSuperEllipse.qMean;           
            xv = linspace(xmin,xmax,n);
            yv = linspace(ymin,ymax,n);
            [xq,yq] = meshgrid(xv,yv); 
            vq = griddata(x,y,v,xq,yq);            
        end
        
        function plotQmeanMxMy(obj)
            x = obj.optimalSuperEllipse.mxV;
            y = obj.optimalSuperEllipse.myV;
            v = obj.optimalSuperEllipse.qMean;    
            f = figure();
            n = 50; 
            xmin = 0.011;0.009;
            xmax = 0.989;0.991;
            ymin = 0.011;0.009;
            ymax = 0.989;0.991;           
            [xq,yq,vq] = obj.interpolateQmeanMxMy(n,xmin,xmax,ymin,ymax);
            h = surf(xq,yq,vq);
            hold on
            plot3(x,y,v,'+')
            xlabel('$m_1$','Interpreter','latex');
            ylabel('$m_2$','Interpreter','latex');
            zlabel('q');
            tN = '\textrm{Optimal SuperEllipse exponent}';
            title(['$',tN,'$'],'interpreter','latex')
            fp = surfPrinter(f);            
            filePath = [obj.outputPath,'QmeanMxMy'];
            fp.print(filePath);
        end
        
    end
    

end