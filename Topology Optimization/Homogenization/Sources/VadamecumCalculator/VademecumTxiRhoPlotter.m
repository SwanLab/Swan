classdef VademecumTxiRhoPlotter < VademecumPlotter
    
    properties (SetAccess = private, GetAccess = public)        
        feasibleIndex
    end
    
    properties (Access = protected)
        XYname = ' TxiRho'
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
            ind = obj.feasibleIndex;
            obj.xV = obj.chi(ind);
            obj.yV = obj.volume(ind);            
        end
        
        function plot(obj)
            obj.plotMxMy();
            obj.plotHomogenizedTensor();
            obj.plotAmplificatorTensor();
        end
        
    end
    
    methods (Access = private)
        
        function computeTxiMxMyVariables(obj)
            for i = 1:length(obj.mxV)
                for j = 1:length(obj.myV)
                    mx = obj.mxV(i);
                    my = obj.myV(j);
                    obj.chi(i,j) = mx/my;
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
            obj.printWidthVariable(obj.mxT,'Mx');
            obj.printWidthVariable(obj.myT,'My');
        end
        
        function printWidthVariable(obj,val,name)
            obj.fileName    = name;
            obj.value2print = val;
            obj.plotFigure();
            obj.printFigure();
        end
        
    end
    
    methods (Access = protected)
        
        function plotFigure(obj)
            x = obj.xV;
            y = obj.yV;
            z = obj.value2print(obj.feasibleIndex);
            ncolors = 50;
            tri = delaunay(x,y);
            obj.fig = figure;
            tricontour(tri,x,y,z,ncolors)
            colorbar
            hold on
            plot(x,y,'+');
            xlabel('$\frac{m1}{m2}$','Interpreter','latex');
            ylabel('\rho');
            obj.addTitle();
        end
        
    end
    
    
    
end