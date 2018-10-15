classdef Monitoring_LevelSet_3D < Monitoring_LevelSet
    properties
        filter
        rotation_per_it
    end
    
    methods
        function obj = Monitoring_LevelSet_3D(settings,mesh,monitoring_ON,plotting_ON)
            obj@Monitoring_LevelSet(settings,mesh,monitoring_ON,plotting_ON);
            obj.rotation_per_it = settings.rotation_per_it;
            obj.filter = Filter_Boundary.create(settings);
            obj.filter.preProcess;
            obj.filter.computeSurroundingFacets;
        end
        
        function setPlottingFigure(obj)
            figure_position = obj.getPlotFigurePosition;
            set(obj.plotting_figure,'Pointer','arrow','Color',[1 1 1],'Name','Finite Element Model','NumberTitle','off','Position',figure_position);
            axis equal; axis off; hold on;
            fac = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
            lx=max(obj.mesh.coord(:,1));
            ly=max(obj.mesh.coord(:,2));
            lz=max(obj.mesh.coord(:,3));
            patch(axes(obj.plotting_figure),'Faces',fac,'Vertices',[0 0 0; 0 ly 0; lx ly 0; lx 0 0; 0 0 lz; 0 ly lz; lx ly lz; lx 0 lz],'FaceColor','w','FaceAlpha',0.0);
            
            view(30,30);
            set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            axis equal off
        end
        
        function refreshPlottingFigure(obj)
            axis equal; axis off;  hold on;
            fac = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
            lx=max(obj.mesh.coord(:,1));
            ly=max(obj.mesh.coord(:,2));
            lz=max(obj.mesh.coord(:,3));
            patch(axes(obj.plotting_figure),'Faces',fac,'Vertices',[0 0 0; 0 ly 0; lx ly 0; lx 0 0; 0 0 lz; 0 ly lz; lx ly lz; lx 0 lz],'FaceColor','w','FaceAlpha',0.0);
            set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            axis equal off
        end
        
        function plotX(obj,x)
%             iso = 0;
%             load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));
%             for n = 1:length(x)
%                 c(b1(n,1),b1(n,2),b1(n,3)) = x(n);
%             end
%             c=permute(c,[3 2 1]);
%             c(c==0)=-eps;
%             
%             cc=(iso+1000)*ones(size(c)+2);
%             cc(2:end-1,2:end-1,2:end-1)=c;
%             
%             [Y,X,Z]=meshgrid(-dim(1,2)/div(1,2):dim(1,2)/div(1,2):dim(1,2)+dim(1,2)/div(1,2),...
%                 -dim(1,1)/div(1,1):dim(1,1)/div(1,1):dim(1,1)+dim(1,1)/div(1,1),...
%                 -dim(1,3)/div(1,3):dim(1,3)/div(1,3):dim(1,3)+dim(1,3)/div(1,3));
%             
%             [F,V,col] = MarchingCubes(X,Y,Z,cc,iso);

            [boundary_facets_coordinates,boundary_facets_connectivities] = obj.filter.computeBoundaryFacets(x);
            
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

