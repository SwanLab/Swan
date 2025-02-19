classdef StressAndMaxNormWithExponentPlotter < StressNormPlotter
    
    methods (Access = public)
        
        function obj = StressAndMaxNormWithExponentPlotter(d)
            obj.init(d);
        end
        
        function plot(obj)
            obj.fig = figure(1);
            obj.plotSmoothInfo();
            obj.plotNonSmoothInfo();
            obj.addLegend();
            obj.print();
        end
        
    end
    
    methods (Access = private)
        
        function plotSmoothInfo(obj)
            dB = obj.smoothDB;
            obj.plotStressNorm(dB);
            obj.plotStressMax(dB);
        end
        
        function plotNonSmoothInfo(obj)
            dB = obj.nonSmoothDB;
            obj.plotStressNorm(dB);
            obj.plotStressMax(dB);
        end
        
        function plotStressNorm(obj,dB)
            obj.x        = dB.pNorm;
            obj.y        = dB.invStressNorm;
            obj.lineInfo = [dB.color,'-+'];
            obj.appendPlot()
        end
        
        function plotStressMax(obj,dB)
            obj.x            = dB.pNorm;
            obj.y            = dB.invStressMax*ones(1,length(obj.x));
            obj.lineInfo     = [dB.color,'--'];
            obj.appendPlot()
        end
        
        function addLegend(obj)
            legStr = {'Smooth','NonSmooth','SmoothMaxNorm','NonSmoothMaxNorm'};
            lines = obj.plotLines;
            nlinesType = length(legStr);
            ilines = 1:obj.nExp:obj.nExp*nlinesType;
            legLines = [lines{ilines(1)} lines{ilines(2)} lines{ilines(3)} lines{ilines(4)}];
            legend(legLines,legStr);
        end
        
        function print(obj)
            sFig = obj.fig;
            ca = get(sFig,'CurrentAxes');
            xl = get(ca,'xlabel');
            set(xl,'string','Pnorm');            
            figureName = 'StressNormWithPnorm';
            outPutFigName = fullfile(obj.outPutPath,figureName);
            pL   = obj.plotLines;
            fp   = plotPrinter(sFig,pL);
            fp.print(outPutFigName)
        end
        
    end
    
end