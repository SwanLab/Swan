classdef BoundayConditionsPlotter_3D < BoundayConditionsPlotter_Abstract
    
    methods (Access = public)
        
        function obj = BoundayConditionsPlotter_3D(axes,mesh,bc)
            obj.init(axes,mesh,bc);            
        end        
        
        function plotDirichlet(obj)
            plot3(obj.axes,obj.mesh.coord(obj.iD,1),obj.mesh.coord(obj.iD,2),obj.mesh.coord(obj.iD,3),'>','Color',obj.colorD,'MarkerSize',4,'MarkerFaceColor',obj.colorD)
            quiver3(obj.axes,obj.mesh.coord(obj.iD,1),obj.mesh.coord(obj.iD,2),obj.mesh.coord(obj.iD,3),obj.vD(:,1),obj.vD(:,2),obj.vD(:,3),'Color',obj.colorD,'AutoScaleFactor',obj.scaleD,'LineWidth',obj.lineWidth,'MaxHeadSize',obj.maxHeadSize);
        end
        
        function plotNeumann(obj)
            plot3(obj.axes,obj.mesh.coord(obj.iN,1),obj.mesh.coord(obj.iN,2),obj.mesh.coord(obj.iN,3),'o','MarkerFaceColor',obj.colorN,'Color',obj.colorN)
            quiver3(obj.axes,obj.mesh.coord(obj.iN,1),obj.mesh.coord(obj.iN,2),obj.mesh.coord(obj.iN,3),obj.vN(:,1),obj.vN(:,2),obj.vN(:,3),'Color',obj.colorN,'AutoScaleFactor',obj.scaleN,'LineWidth',obj.lineWidth,'MaxHeadSize',obj.maxHeadSize);
        end
        
    end
    
end