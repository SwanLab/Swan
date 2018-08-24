classdef Monitoring_2D < Monitoring
    methods
        function obj = Monitoring_2D(settings,mesh,monitoring_ON,plotting_ON)
            obj@Monitoring(settings,mesh,monitoring_ON,plotting_ON);
            obj.ndim = 2;
        end
        
        function setPlottingFigure(obj)
            mp = get(0, 'MonitorPositions');
            select_screen = 1;
            if size(mp,1) < select_screen
                select_screen = size(mp,1);
            end
            width = mp(1,3);
            height = mp(1,4);
            size_screen_offset = round([0.7*width,0.52*height,-0.71*width,-0.611*height],0);
            set(obj.plotting_figure,'Position',mp(select_screen,:) + size_screen_offset);
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

