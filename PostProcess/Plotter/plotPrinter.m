classdef plotPrinter < figurePrinter
    
    properties (Access = private)
        %lineWidth = 2;
        lineWidth = 4;
        pObj
    end
    
    methods (Access = public)
        
        function obj = plotPrinter(figureID,plotObject)
            obj.init(figureID,plotObject);
            obj.setFontSizes();
            obj.setPlotsLineWidth()
            obj.setPositionAndSize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,figureID,plotObject)
            obj.figID = figureID;
            obj.pObj = plotObject;
        end
        
        function setPlotsLineWidth(obj)
            for iplots = 1:numel(obj.pObj)
                set(obj.pObj{iplots},'LineWidth',obj.lineWidth);
            end
        end
        
    end
    
    
end