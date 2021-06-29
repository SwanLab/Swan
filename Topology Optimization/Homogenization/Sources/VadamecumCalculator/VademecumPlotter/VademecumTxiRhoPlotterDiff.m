classdef VademecumTxiRhoPlotterDiff < VademecumPlotter
    
    properties (SetAccess = private, GetAccess = public)
        feasibleIndex
    end
    
    properties (Access = protected)
       name = 'TxiRho'; 
    end    
        
    properties (Access = private)
        chi
        smoothDB        
        nonSmoothDB
        txiN
        rhoN
        iN
        txiS
        rhoS
        iS
        txiV
        rhoMin
        rhoMax
        indexSort
        plotter
    end
    
    methods (Access = public)
        
        function obj = VademecumTxiRhoPlotterDiff(d)
            obj.init(d);
            obj.computeTxiVariable();
            obj.txiN = obj.chi(obj.iN);
            obj.rhoN = obj.nonSmoothDB.volume(obj.iN);
            obj.txiS = obj.chi(obj.iS);
            obj.rhoS = obj.smoothDB.volume(obj.iS);            
        end
        
        function plot(obj)
            obj.plotHomogenizedTensor();
            obj.plotAmplificatorTensor();
        end
        
    end
    
    methods (Access = private)
        
        function computeTxiVariable(obj)
            for i = 1:length(obj.mxV)
                for j = 1:length(obj.myV)
                    mx = obj.mxV(i);
                    my = obj.myV(j);
                    obj.chi(i,j) = atan(mx/my);
                end
            end
        end
              
        function zDif = computeDiff(obj,varN,varS)
            zN = varN(obj.iN);            
            zS = varS(obj.iS);
            x_c(:,1) = obj.txiS;
            y_c(:,1) = obj.rhoS;
            z_c(:,1) = zS;
            F = scatteredInterpolant(x_c,y_c,z_c);
            ZsInt = F(obj.txiN,obj.rhoN);
            zDif = (ZsInt - zN);
        end                        
        
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.index = [1 1; 2 2;3 3; 1 2;2 3;1 3];            
            obj.smoothDB    = d.smoothDB;
            obj.nonSmoothDB = d.nonSmoothDB;
            obj.mxV = obj.nonSmoothDB.mxV;
            obj.myV = obj.nonSmoothDB.myV;
            obj.iS = intersect(d.iS,d.iN);
            obj.iN = intersect(d.iS,d.iN);
            obj.hasToPrint = d.hasToPrint;
            obj.outPutPath = d.outPutPath;
            obj.microName  = d.microName;            
        end        
        
        function plotHomogenizedTensor(obj)
            obj.tensor{1} = obj.smoothDB.C;
            obj.tensor{2} = obj.nonSmoothDB.C;            
            obj.tensorCase = 'C_';
            obj.tensorTitleName = '\Delta C_';            
            obj.plotTensor();
        end
                
        function plotAmplificatorTensor(obj)
            obj.tensor{1} = obj.smoothDB.invP;
            obj.tensor{2} = obj.nonSmoothDB.invP;  
            obj.tensorCase = 'P_';            
            obj.tensorTitleName = '(\Delta  P^{-1})_';
            obj.plotTensor();
        end     
        
        function obtainTensorComponent(obj)
            i = obj.iIndex;
            j = obj.jIndex;
            tS = obj.tensor{1}(i,j,:,:);
            tS = squeeze(tS);
            tN = obj.tensor{2}(i,j,:,:);
            tN = squeeze(tN);
            Tij = obj.computeDiff(tN,tS);
            obj.value2print = Tij;
        end        
        
        function plotFigure(obj)
            obj.computeRhoMinRhoMaxAndTxiV();            
            obj.plotDifference();
%             obj.plotRhoMinRhoMaxLines();            
            %obj.plotSmoothSuperEllipsePoints();
%           obj.plotRectanglePoints();
            %obj.addTitleAndAxisNames();
            obj.addZeroLevelSetValue();
        end        
        
        function computeRhoMinRhoMaxAndTxiV(obj)
            [obj.txiV,obj.indexSort] = sort(obj.txiS);
            s.q = 2;
            s.txi = obj.txiV;
            s.mxMin = min(obj.smoothDB.mxV);
            s.myMin = min(obj.smoothDB.myV);  
            s.mxMax = max(obj.smoothDB.mxV);            
            s.myMax = max(obj.smoothDB.myV);    
            rhoBounds = SuperEllipseRhoBoundsComputer(s);
            [rMin,rMax] = rhoBounds.compute();      
            obj.rhoMin = rMin;
            obj.rhoMax = rMax;
        end        
        
        function plotDifference(obj)
            x(:,1) = obj.txiS;
            y(:,1) = obj.rhoS;            
            z(:,1) = obj.value2print;
%             ind = ~isnan(z);
%             ncolors = 50;
%             tri = delaunay(x(ind),y(ind));
%             obj.fig = figure;            
%             tricontour(tri,x(ind),y(ind),z(ind),ncolors) 
%             colorbar
                 
            s.fileName = fullfile(obj.outPutPath,[obj.fileName,'XiRho']);
            s.title    = obj.titleName;
            s.axisAdder = XiRhoAxisAdder();
            p  = SuperEllipseExponentContourPlotter(s); 
            p.plot(x,y,z);
%           
            obj.plotter = p;
            
            hold on            
        end
        
        function plotSmoothSuperEllipsePoints(obj)
            x = obj.txiS;
            y = obj.rhoS;             
            plot(x,y,'+');
        end
        
        function plotRectanglePoints(obj)
           rho = obj.rhoN(obj.indexSort);
           txi = obj.txiN(obj.indexSort);
           isLowerRhoMax = rho <= obj.rhoMax;
           isUpperRhoMin = rho >= obj.rhoMin;
           isFeasible = isLowerRhoMax & isUpperRhoMin;
           plot(txi(isFeasible),rho(isFeasible),'+')            
        end
        
        function addTitleAndAxisNames(obj)
            xlabel('$\xi$','Interpreter','latex');
            ylabel('\rho');
            tN = obj.titleName;
            title(['$',tN,'$'],'interpreter','latex')     
            ylim([0 1])   
            set(gca,'xtick',[0:pi/8:pi/2]) 
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})            
        end
        
        function addZeroLevelSetValue(obj)
            x(:,1) = obj.txiS;
            y(:,1) = obj.rhoS;            
            z(:,1) = obj.value2print;        
%             ind = ~isnan(z);
%             ind(:) = true;
%             connec = delaunay(x(ind),y(ind));     
%             s.x = x;
%             s.y = y;
%             s.connec = connec;
%             c = ConnecWithQualityComputer(s);
%             connec = c.compute();                        
%             v = [0,0];
%             [~,c] = tricontour(connec,x(ind),y(ind),z(ind),v);
%             
            
%             
%             
%             [xi, yi] = meshgrid(...
%             linspace(min(x),max(x)),...
%             linspace(min(y),max(y)));
%             zi = griddata(x,y,z, xi,yi);
            x2 = reshape(x,20,20);
            y2 = reshape(y,20,20);
            z2 = reshape(z,20,20);

           [M,c] = contourf(x2,y2,z2,[0 0]);
           c.Fill = false;
           c.LineWidth = 3;
           c.EdgeColor = 'r';  
           obj.plotter.print();
        end
        
        function plotRhoMinRhoMaxLines(obj)
           obj.plotRhoLine(obj.txiV,obj.rhoMin);
           obj.plotRhoLine(obj.txiV,obj.rhoMax);
        end
        
        function plotRhoLine(obj,txiV,rho)
            h = plot(txiV,rho,['-','k']);
            set(h,'LineWidth',2);        
        end
        
    end
    
    
    
end