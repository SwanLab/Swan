classdef VademecumPlotter < handle
    
    properties (Access = private)
        mxV
        myV
        C
        invP
        index
        iIndex
        jIndex
        tensor
        tensorComp
        tensorName
        tensorCase
        microName
        outPutPath        
        hasToPrint
        fig
    end
    
    
    methods (Access = public)
        
        function obj = VademecumPlotter(d)
            obj.init(d);
        end
        
        function plot(obj)
            obj.plotHomogenizedTensor();
            obj.plotAmplificatorTensor();            
        end
    end
    
    
    methods (Access = private)
        
        function init(obj,d)
            obj.index = [1 1; 2 2;3 3; 1 2;2 3;1 3];
            obj.mxV        = d.mxV;
            obj.myV        = d.myV;
            obj.C          = d.C;
            obj.invP       = d.invP;
            obj.hasToPrint = d.hasToPrint;
            obj.outPutPath = d.outPutPath;
            obj.microName  = d.microName;
        end
        

        
        function plotHomogenizedTensor(obj)
            obj.tensor = obj.C;
            obj.tensorCase = 'C_';
            obj.plotTensor();
        end
        
        function plotAmplificatorTensor(obj)
            obj.tensor = obj.invP;
            obj.tensorCase = 'Pinv_';
            obj.plotTensor();
        end
        
        function plotTensor(obj)
            nPlot = length(obj.index);
            for iplot = 1:nPlot
                obj.obtainIJindex(iplot);
                obj.obtainTensorComponent();
                obj.createTensorName();
                obj.plotFigure();
                obj.printFigure();
            end
        end
        
        function obtainIJindex(obj,iplot)
            obj.iIndex = obj.index(iplot,1);
            obj.jIndex = obj.index(iplot,2);
        end
        
        function obtainTensorComponent(obj)
            i = obj.iIndex;
            j = obj.jIndex;
            Tij = obj.tensor(i,j,:,:);
            Tij = squeeze(Tij);
            obj.tensorComp = Tij;
        end
        
        function createTensorName(obj)
            i = obj.iIndex;
            j = obj.jIndex;            
            obj.tensorName = [obj.tensorCase,'{',num2str(i),num2str(j),'}'];
        end
        
        function addTitle(obj)
            tN = obj.tensorName;
            title([tN,' for ',obj.microName])
        end
        
        function plotFigure(obj)
            obj.fig = figure();
            x = obj.mxV;
            y = obj.myV;
            z = obj.tensorComp;
            contour(x,y,z,50);
            xlabel('mx');
            ylabel('my');
            obj.addTitle();
            colorbar;
        end
        
        function printFigure(obj)
            if obj.hasToPrint
                outPutFigName = [obj.outPutPath,obj.tensorName];
                fp = contourPrinter(obj.fig);
                fp.print(outPutFigName)                
            end
        end
        
    end
    
end