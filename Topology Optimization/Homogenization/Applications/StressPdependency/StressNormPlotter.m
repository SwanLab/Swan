classdef StressNormPlotter < handle
    
    properties (Access = protected)
        outPutPath 
        smoothDB
        nonSmoothDB        
        nExp
        indexLines
        x
        y
        lineInfo
        fig
        plotLines        
    end
    
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.smoothDB    = d.smoothDB;
            obj.nonSmoothDB = d.nonSmoothDB;
            obj.nExp        = d.nExp;
            obj.outPutPath  = d.outPath;
            obj.indexLines = 0;
            obj.creatingColorInfo();
        end
        
        function creatingColorInfo(obj)
            obj.smoothDB.color    = 'b';
            obj.nonSmoothDB.color = 'r';
        end        
        
        function appendPlot(obj)
            figure(obj.fig);
            hold on
            p = plot(obj.x,obj.y,obj.lineInfo);
            obj.savePlotLines(p);
        end
        
        function savePlotLines(obj,p)
            lines = obj.obtainLines();
            for iline = 1:length(p)
                obj.plotLines{lines(iline)} = p(iline);
            end
        end
        
        function lines = obtainLines(obj)
            nLines = size(obj.y,1);
            first = obj.indexLines + 1;
            last  = obj.indexLines + nLines;
            lines = first:last;
            obj.indexLines = last;
        end        
        
    end
    
end
