classdef Monitoring_Else < Monitoring
    methods
        function obj = Monitoring_Else(settings,mesh,monitoring_ON,plotting_ON)
            obj@Monitoring(settings,mesh,monitoring_ON,plotting_ON);
            obj.ndim = 2;
        end
        
        function setPlottingFigure(obj)
            figure_position = obj.getPlotFigurePosition;
            set(obj.plotting_figure,'Pointer','arrow','Color',[1 1 1],'Name','Finite Element Model','NumberTitle','off','Position',figure_position);
            %                 obj.plotting_figure = trisurf(obj.mesh.connec,obj.mesh.coord(:,1),obj.mesh.coord(:,2),obj.mesh.coord(:,3),double(rho_nodal), ...
            %                     'EdgeColor','none','LineStyle','none','FaceLighting','phong');
            
            obj.plotting_figure=patch('Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,'FaceVertexCData',zeros(size(obj.mesh.coord,1),1),'FaceColor','flat',...
                'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
            
            colormap(flipud(gray));
            set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            
            axis equal
            axis off
        end
        
        function plotX(obj,x)
            if any(x<0)
                rho_nodal = x<0;
            else
                rho_nodal = x;
            end
            
            set(obj.plotting_figure,'FaceVertexCData',double(rho_nodal));
%             set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            
            if obj.showBC
                obj.plotBoundaryConditions;
            end
            drawnow;
        end
    end
end

