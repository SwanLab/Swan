classdef figurePrinter < handle
    
    properties (Access = protected)
        x0 = 1920;
        y0 = 500;
        width = 600;
        height = 600;
        %labelFontSize = 30;
        labelFontSize = 200;
        %fontSize = 20;
        fontSize = 30;        
        outFormat = '-dpng';
        figID
    end
   
    methods (Access = public)
        
        function print(obj,figName)
            obj.printFigure(figName)
        end
        
    end    
    
    methods (Access = protected)        
        
        function setFontSizes(obj)
            obj.setCurrentAxisFontSize();
            obj.setXlabelFontSize();
            obj.setYlabelFontSize();
        end        
        
        function setPositionAndSize(obj)
            set(obj.figID,'units','points','position',[obj.x0,obj.y0,obj.width,obj.height])
        end
        
        function printFigure(obj,figName)
            print(obj.figID,figName,obj.outFormat);
        end
        
    end
    
    methods (Access = private)
        
        function setCurrentAxisFontSize(obj)
            set(gca,'fontsize',obj.fontSize)
        end
        
        function setXlabelFontSize(obj)
            ca = get(obj.figID,'CurrentAxes');
            xl = get(ca,'xlabel');
            set(xl,'FontSize',obj.labelFontSize)
        end
        
        function setYlabelFontSize(obj)
            ca = get(obj.figID,'CurrentAxes');
            yl = get(ca,'ylabel');
            set(yl,'FontSize',obj.labelFontSize)
        end                
        
    end
    
    
end