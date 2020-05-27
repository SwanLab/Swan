classdef BoundayConditionsPlotter_2D < BoundayConditionsPlotter_Abstract
    
    methods (Access = public)
        
        function obj = BoundayConditionsPlotter_2D(axes,mesh,bc)
            obj.init(axes,mesh,bc);            
        end
        
        function plotDirichlet(obj)
            if ~isempty(obj.iD)
                plot(obj.axes,obj.mesh.coord(obj.iD,1),obj.mesh.coord(obj.iD,2),'>','Color',obj.colorD,'MarkerSize',4,'MarkerFaceColor',obj.colorD)
                quiver(obj.axes,obj.mesh.coord(obj.iD,1),obj.mesh.coord(obj.iD,2),obj.vD(:,1),obj.vD(:,2),'Color',obj.colorD,'AutoScaleFactor',obj.scaleD,'LineWidth',obj.lineWidth,'MaxHeadSize',obj.maxHeadSize);
            end
        end
        
        function plotNeumann(obj)
            if ~isempty(obj.iN)
                plot(obj.axes,obj.mesh.coord(obj.iN,1),obj.mesh.coord(obj.iN,2),'o','MarkerFaceColor',obj.colorN,'Color',obj.colorN)
                quiver(obj.axes,obj.mesh.coord(obj.iN,1),obj.mesh.coord(obj.iN,2),obj.vN(:,1),obj.vN(:,2),'Color',obj.colorN,'AutoScaleFactor',obj.scaleN,'LineWidth',obj.lineWidth,'MaxHeadSize',obj.maxHeadSize);
            end
        end
        
    end
    
end