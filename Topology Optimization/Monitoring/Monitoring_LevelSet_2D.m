classdef Monitoring_LevelSet_2D < Monitoring_LevelSet
    methods
        function obj = Monitoring_LevelSet_2D(settings,mesh,monitoring_ON,plotting_ON)
            obj@Monitoring_LevelSet(settings,mesh,monitoring_ON,plotting_ON);
        end
        
        function setPlottingFigure(obj)
            figure_position = obj.getPlotFigurePosition;
            set(obj.plotting_figure,'Pointer','arrow','Color',[1 1 1],'Name','Finite Element Model','NumberTitle','off','Position',figure_position);
            axis equal; axis off; hold on;
            set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            axis equal off
        end
        
        function refreshPlottingFigure(obj)
            set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            axis equal off
        end
        
        function plotX(obj,x)
            obj.unfitted_mesh = Mesh_Unfitted_2D(obj.mesh.duplicate,obj.geometry.interpolation,x);
            obj.unfitted_mesh.computeCutMesh;
            
            [boundary_facets_coordinates,boundary_facets_connectivities] = obj.unfitted_mesh.computeBoundaryFacets(x);
            
            set(0, 'CurrentFigure', obj.plotting_figure)
            [az,el] = view;
            clf
            hold on
            
            obj.refreshPlottingFigure;
            view([az+obj.rotation_per_it,el]);
            
            patch('vertices',boundary_facets_coordinates,'faces',boundary_facets_connectivities,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
            light
            axis equal off
            
            if obj.showBC
                obj.plotBoundaryConditions;
            end
            drawnow;
        end
    end
    
    methods (Static)
        function figure_position = getPlotFigurePosition
            figure_position = get(0, 'Screensize');
        end
    end
end

