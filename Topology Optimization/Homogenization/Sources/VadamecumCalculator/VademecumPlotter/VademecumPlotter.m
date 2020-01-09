classdef VademecumPlotter < handle
        
    properties (Access = protected, Abstract)
       name 
    end
    
    properties (Access = protected)
        index  
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
        titleName
        tensorTitleName
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
                       
        function printFigure(obj)
            if obj.hasToPrint
                outPutFigName = [obj.outPutPath,obj.fileName,obj.name];
                fp = contourPrinter(obj.fig);
                fp.print(outPutFigName)
            end
        end   
        
        function plotHomogenizedTensor(obj)
            obj.tensor = obj.C;
            obj.tensorCase = 'C_';
            obj.tensorTitleName = 'C_';
            obj.plotTensor();
        end
        
        function plotAmplificatorTensor(obj)
            obj.tensor = obj.invP;
            obj.tensorCase = 'Pinv_';
            obj.tensorTitleName = '(P^{-1})_';
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
        
        function obtainTensorComponent(obj)
            i = obj.iIndex;
            j = obj.jIndex;
            Tij = obj.tensor(i,j,:,:);
            Tij = squeeze(Tij);
            obj.value2print = Tij;
        end       
        
        function plotHomogenizedTensorIsotropy(obj)
            obj.fileName = 'Cisotropy';
            obj.titleName = '\textrm{Isotropy of C}';
            obj.value2print = obj.computeIsotropyNorm(obj.C);
            obj.plotFigure();
            obj.printFigure();
        end        
        
    end
    
    methods (Access = private)
        
        function obtainIJindex(obj,iplot)
            obj.iIndex = obj.index(iplot,1);
            obj.jIndex = obj.index(iplot,2);
        end
        
        function createTensorName(obj)
            i = obj.iIndex;
            j = obj.jIndex;
            obj.fileName  = [obj.tensorCase,'{',num2str(i),num2str(j),'}'];
            obj.titleName = [obj.tensorTitleName,'{',num2str(i),num2str(j),'}'];
        end               
        
    end
    
    methods (Access = protected, Abstract)
        plotFigure(obj)
    end
    
     methods (Access = protected, Static)
        
        function in = computeIsotropyNorm(tensor)
            t11 = tensor(1,1,:,:); 
            t22 = tensor(2,2,:,:); 
            t12 = tensor(1,2,:,:);
            t33 = tensor(3,3,:,:);            
            t13 = tensor(1,3,:,:);
            t23 = tensor(2,3,:,:);                        
            eq(1,:,:) = abs(t11 - t22);
            eq(2,:,:) = abs(t33 - 0.5*(t11 + t22) - t12);
            eq(3,:,:) = abs(t23);
            eq(4,:,:) = abs(t13);
            nt = abs(t11) + abs(t22) + abs(t12) + abs(t33) + abs(t13) + abs(t23);
            nt = squeeze(nt);
            in = squeeze(max(eq,[],1));
            in = in./nt;
        end        
        
    end
    
end