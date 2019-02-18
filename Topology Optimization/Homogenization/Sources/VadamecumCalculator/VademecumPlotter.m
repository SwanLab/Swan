classdef VademecumPlotter < handle
    
    properties (Abstract, Access = protected)
        XYname
    end
    
    properties (Access = protected)
        index = [1 1; 2 2;3 3; 1 2;2 3;1 3]; 
        iIndex
        jIndex
        value2print   
        microName
        outPutPath   
        hasToPrint     
        fileName
        fig    
        tensor
        tensorCase   
        xV
        yV   
        mxV
        myV
        C
        invP
        volume 
    end
        
    methods (Access = protected)
        
        function init(obj,d)
            obj.index = [1 1; 2 2;3 3; 1 2;2 3;1 3];
            obj.mxV        = d.mxV;
            obj.myV        = d.myV;
            obj.C          = d.C;
            obj.invP       = d.invP;
            obj.volume     = d.volume;
            obj.hasToPrint = d.hasToPrint;
            obj.outPutPath = d.outPutPath;
            obj.microName  = d.microName;
        end        
        
        function addTitle(obj)
            fN = obj.fileName;
            title([fN,' for ',obj.microName])
        end        
        
        function printFigure(obj)
            if obj.hasToPrint
                outPutFigName = [obj.outPutPath,obj.fileName];
                fp = contourPrinter(obj.fig);
                fp.print(outPutFigName)
            end
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
        
    end
    
    methods (Access = private)
        
        function obtainIJindex(obj,iplot)
            obj.iIndex = obj.index(iplot,1);
            obj.jIndex = obj.index(iplot,2);
        end
        
        function obtainTensorComponent(obj)
            i = obj.iIndex;
            j = obj.jIndex;
            Tij = obj.tensor(i,j,:,:);
            Tij = squeeze(Tij);
            obj.value2print = Tij;
        end
        
        function createTensorName(obj)
            i = obj.iIndex;
            j = obj.jIndex;
            obj.fileName = [obj.tensorCase,'{',num2str(i),num2str(j),'}',obj.XYname];
        end               
        
    end
    
    methods (Access = protected, Abstract)
        plotFigure(obj)
    end
    
end