classdef Monitoring_LevelSet_3D < Monitoring
    properties
        filter
    end
    
    methods
        function obj = Monitoring_LevelSet_3D(settings,mesh,monitoring_ON,plotting_ON)
            obj@Monitoring(settings,mesh,monitoring_ON,plotting_ON);
            obj.ndim = 3;
            obj.filter =  Filter.create(settings);
            obj.filter.preProcess;
            obj.filter.computeSurroundingFacets;
        end
        
        function setPlottingFigure(obj)
            set(obj.plotting_figure,'Pointer','arrow','Color',[1 1 1],'Name','Finite Element Model','NumberTitle','off');
            axis equal; axis off; view(3); hold on;
            fac = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
            lx=max(obj.mesh.coord(:,1));
            ly=max(obj.mesh.coord(:,2));
            lz=max(obj.mesh.coord(:,3));
            patch(axes(obj.plotting_figure),'Faces',fac,'Vertices',[0 0 0; 0 ly 0; lx ly 0; lx 0 0; 0 0 lz; 0 ly lz; lx ly lz; lx 0 lz],'FaceColor','w','FaceAlpha',0.0);
            
            rotate3d(gca);
            set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            view(30,30);
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
            clf
            hold on
            
            obj.setPlottingFigure;
            
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
end

