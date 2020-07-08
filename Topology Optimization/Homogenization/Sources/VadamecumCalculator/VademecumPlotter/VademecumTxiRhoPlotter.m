classdef VademecumTxiRhoPlotter < VademecumPlotter
    
    properties (SetAccess = private, GetAccess = public)        
        feasibleIndex
    end
    
    properties (Access = protected)
       name = 'TxiRho'; 
    end
    
    properties (Access = private)
        mxT
        myT
        chi    
    end
    
    methods (Access = public)
        
        function obj = VademecumTxiRhoPlotter(d)
            obj.init(d);
            obj.computeTxiMxMyVariables();
            obj.computeFeasibleIndex();
            obj.feasibleIndex = 1:length(obj.mxV)*length(obj.myV);
            ind = obj.feasibleIndex;
            obj.xV = obj.chi(ind);
            obj.yV = obj.volume(ind);            
        end
        
        function plot(obj)
            obj.plotMxMy();
            obj.plotHomogenizedTensor();
            obj.plotHomogenizedTensorIsotropy();            
            obj.plotAmplificatorTensor();
        end
        
    end
    
    methods (Access = private)
        
        function computeTxiMxMyVariables(obj)
            for i = 1:length(obj.mxV)
                for j = 1:length(obj.myV)
                    mx = obj.mxV(i);
                    my = obj.myV(j);
                    obj.chi(i,j) = atan(mx/my);
                    obj.mxT(i,j) = mx;
                    obj.myT(i,j) = my;
                end
            end
        end
                  
        function computeFeasibleIndex(obj)
            d.mx = obj.mxV;
            d.my = obj.myV;
            d.chi = obj.chi;
            d.rho = obj.volume;
            fC = FeasibleIndexComputer(d);
            obj.feasibleIndex = fC.index;
        end
        
        function plotMxMy(obj)
            obj.printWidthVariable(obj.mxT,'m_1');
            obj.printWidthVariable(obj.myT,'m_2');
        end
        
        function printWidthVariable(obj,val,name)
            obj.fileName    = name;
            obj.titleName   = name;
            obj.value2print = val;
            obj.plotFigure();
            obj.printFigure();
        end
        
    end
    
    methods (Access = protected)
        
        function plotFigure(obj)
            x(:,1) = obj.xV;
            y(:,1) = obj.yV;
            z(:,1) = obj.value2print(obj.feasibleIndex);
            
            
            s.fileName = fullfile(obj.outPutPath,[obj.fileName,'XiRho']);
            s.title    = obj.titleName;
            s.axisAdder = XiRhoAxisAdder();
            p  = SuperEllipseExponentContourPlotter(s); 
            p.plot(x,y,z);
%             
%             ncolors = 50;
%             tri = delaunay(x,y);
%             obj.fig = figure;
%             tricontour(tri,x,y,z,ncolors)
%             colorbar
%             hold on
%             plot(x,y,'+');
%             ylim([0 1])
%             xlabel('$\xi$','Interpreter','latex');
%             ylabel('\rho');
%             tN = obj.titleName;
%             title(['$',tN,'$'],'interpreter','latex')
%             set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
%             set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})            
        end
        
    end
    
    
    
end