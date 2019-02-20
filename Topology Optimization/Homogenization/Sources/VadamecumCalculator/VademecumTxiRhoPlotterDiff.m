classdef VademecumTxiRhoPlotterDiff < VademecumPlotter
    
    properties (SetAccess = private, GetAccess = public)
        feasibleIndex
    end
    
    properties (Access = protected)
        XYname = ' TxiRho'
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
                    obj.chi(i,j) = mx/my;
                end
            end
        end
              
        function zDif = computeDiff(obj,varN,varS)
            zN = varN(obj.iN);            
            zS = varS(obj.iS);
            ZnInt = griddata(obj.txiN,obj.rhoN,zN,obj.txiS,obj.rhoS);            
            zDif = abs(zS - ZnInt);            
        end                        
        
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.smoothDB    = d.smoothDB;
            obj.nonSmoothDB = d.nonSmoothDB;
            obj.mxV = obj.nonSmoothDB.mxV;
            obj.myV = obj.nonSmoothDB.myV;
            obj.iS = d.iS;
            obj.iN = d.iN;
            obj.hasToPrint = d.hasToPrint;
            obj.outPutPath = d.outPutPath;
            obj.microName  = d.microName;            
        end        
        
        function plotHomogenizedTensor(obj)
            obj.tensor{1} = obj.smoothDB.Ctensor;
            obj.tensor{2} = obj.nonSmoothDB.Ctensor;            
            obj.tensorCase = 'C_';
            obj.plotTensor();
        end
                
        function plotAmplificatorTensor(obj)
            obj.tensor{1} = obj.smoothDB.PinvTensor;
            obj.tensor{2} = obj.nonSmoothDB.PinvTensor;                
            obj.tensorCase = 'Pinv_';
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
            x = obj.txiN;
            y = obj.rhoN;
            z = obj.value2print;
            ind = ~isnan(z);
            ncolors = 50;
            tri = delaunay(x(ind),y(ind));
            obj.fig = figure;            
            tricontour(tri,x(ind),y(ind),z(ind),ncolors) 
            colorbar
            hold on
            plot(x,y,'+');
            plot(obj.txiS,obj.rhoS,'+')            
            xlabel('$\frac{m1}{m2}$','Interpreter','latex');
            ylabel('\rho');
            obj.addTitle();
        end        
        
    end
    
    
    
end