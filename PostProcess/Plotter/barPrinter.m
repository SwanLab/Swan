classdef barPrinter < figurePrinter
    
    properties (Access = private)
        pObj
    end
    
    methods (Access = public)
        
        function obj = barPrinter(figureID,plotObject)
            obj.init(figureID,plotObject);
            obj.setFontSizes();
            obj.setPositionAndSize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,figureID,plotObject)
            obj.figID = figureID;
            obj.pObj = plotObject;
            obj.labelFontSize = 10;
            obj.fontSize = 10;
        end
        
        function setPlotsLineWidth(obj)
            for iplots = 1:numel(obj.pObj)
                set(obj.pObj{iplots},'LineWidth',obj.lineWidth);
            end
        end
        
    end
    
    
end