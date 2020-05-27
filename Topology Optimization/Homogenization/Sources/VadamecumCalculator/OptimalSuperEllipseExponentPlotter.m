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
            obj.plotQdesvTxiRho();
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
        
        function plotQdesvTxiRho(obj)
            obj.plotQdesvTxiRhoContour();  
            obj.plotQdesvTxiRhoTrisurf();                                    
        end
        
        function plotContour(obj,q,fName)
            
            x = obj.optimalSuperEllipse.txiV;
            y = obj.optimalSuperEllipse.rhoV;
            z = q;
            tri = delaunay(x,y);
            
            s.coord = [x',y'];
            s.connec = tri;
            m = Mesh().create(s);
            qua = m.computeElementQuality';
            isQ = qua > 0.02;
            f = figure();            
            h = trisurf(tri(isQ,:),x,y,z);
            colorbar
            h.FaceColor = 'interp';
            %h.EdgeColor = 'none';
            view(2)
            
            ylim([0 1])
            xlabel('$\xi$','Interpreter','latex');
            ylabel('\rho');
            % tN = '\textrm{Optimal SuperEllipse exponent}';
            % title(['$',tN,'$'],'interpreter','latex')
            set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})
            fp = contourPrinter(f);
            filePath = [obj.outputPath,fName];
            fp.print(filePath);
                        
        end
        
        function plotQmeanTxiRhoContour(obj)
            q = obj.optimalSuperEllipse.qMean;            
            fName = 'QmeanTxiRhoContour';
            obj.plotContour(q,fName);       
        end
        
        function plotQdesvTxiRhoContour(obj)
            q = obj.optimalSuperEllipse.qDesv;            
            fName = 'QdesvTxiRhoContour';
            obj.plotContour(q,fName);                         
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
            q = obj.optimalSuperEllipse.qMean;
            fName = 'QmeanTxiRhoSurf';
            obj.plotTrisurf(q,fName);
        end
        
        function plotQdesvTxiRhoTrisurf(obj)
            q = obj.optimalSuperEllipse.qDesv;
            fName = 'QdesvTxiRhoSurf';
            obj.plotTrisurf(q,fName);            
        end
        
        function plotTrisurf(obj,q,fName)
            x = obj.optimalSuperEllipse.txiV;
            y = obj.optimalSuperEllipse.rhoV;
            z = q; 
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
            filePath = [obj.outputPath,fName];
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