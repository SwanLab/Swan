classdef figurePlotter < handle
    
    properties (Access = private)
        x0 = 1920;
        y0 = 500;
        width = 600;
        height = 600;
        lineWidth = 1.5;
        xlabelFontSize = 30;
        fontSize = 20;
        outFormat = '-dpng';
        figID
        pObj
        xlabelName
    end
    
    methods (Access = public)
        
        function obj = figurePlotter(figureID,plotObject,xlabelName)
            obj.init(figureID,plotObject,xlabelName);
            obj.setCurrentAxisFontSize();
            obj.setXlabelFontSize();
            obj.setPlotsLineWidth()
            obj.setPositionAndSize();
        end
        
        function print(obj,figName)
            obj.printFigure(figName)
        end
        
        
    end
        
    methods (Access = private)
        
        function init(obj,figureID,plotObject,xName)
            obj.figID = figureID;
            obj.pObj = plotObject;
            obj.xlabelName = xName;
        end
        
        function setCurrentAxisFontSize(obj)
            set(gca,'fontsize',obj.fontSize)
        end
        
        function setXlabelFontSize(obj)
            ca = get(obj.figID,'CurrentAxes');
            xl = get(ca,'xlabel');
            set(xl,'string',obj.xlabelName,'FontSize',obj.xlabelFontSize)
        end
        
        function setPlotsLineWidth(obj)
            for iplots = 1:numel(obj.pObj)
                set(obj.pObj{iplots},'LineWidth',1.5);
            end
        end
        
        function setPositionAndSize(obj)
            set(obj.figID,'units','points','position',[obj.x0,obj.y0,obj.width,obj.height])
        end
        
        function printFigure(obj,figName)
            print(obj.figID,figName,obj.outFormat);
        end
    end
    
    
end