classdef StressMaxWithMeshSizePlotter < StressNormPlotter
    
    methods (Access = public)
        
        function obj = StressMaxWithMeshSizePlotter(d)
            obj.init(d);            
        end
        
        function plot(obj)
            obj.fig = figure(2);
            obj.plotSmoothInfo();
            obj.plotNonSmoothInfo();
            obj.addLegend();
            obj.print();            
        end
        
    end
    
    methods (Access = private)
        
        function plotSmoothInfo(obj)
            d = obj.smoothDB;
            obj.plotInfo(d);
        end
        
        function plotNonSmoothInfo(obj)
            d = obj.nonSmoothDB;
            obj.plotInfo(d);
        end
        
        function plotInfo(obj,dB)
            obj.x        = dB.meshSize;
            obj.y        = (1./dB.invStressMax)';
            obj.lineInfo = [dB.color,'-+'];
            obj.appendPlot();
        end
        
        function addLegend(obj)
            legStr = {'Smooth','NonSmooth'};
            lines = obj.plotLines;
            legLines = [lines{1} lines{2}];
            legend(legLines,legStr);
        end
        
        function print(obj)
            sFig = obj.fig;
            ca = get(sFig,'CurrentAxes');
            xl = get(ca,'xlabel');
            set(xl,'string','MeshSize');            
            figureName = 'StressNormWithMeshSize';
            outPutFigName = [obj.outPutPath,figureName];
            pL   = obj.plotLines;
            fp = plotPrinter(sFig,pL,xlabelName);
            fp.print(outPutFigName)
        end
    end
        
end