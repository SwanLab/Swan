classdef DesignVarMonitor_AreaColumn < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarName = 'Area Column'
    end
    
    properties (Access = private)
        polygon
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_AreaColumn(cParams)
            obj@DesignVarMonitor_Abstract(cParams);
        end
        
        function plot(obj)
            scale = 0.3;
            obj.createPolygon(scale);
            obj.plotFigure();
        end
        
    end

    methods (Access = private)

        function createPolygon(obj,scl)
            z = obj.designVar.getColumnArea;
            coord = obj.mesh.coord;
            nnod     = obj.mesh.nelem+1;
            dimFig   = 2;
            vertElem = 4;
            vertex = zeros(vertElem*obj.mesh.nelem+1,dimFig);
            for iNod = 1:nnod-1
                vertex(2*iNod-1,:)   = [scl*z(iNod)/2 coord(iNod)];
                vertex(2*iNod,:) = [scl*z(iNod)/2 coord(iNod+1)];
            end
            vertex = obj.flip(vertex,vertElem,dimFig);
            vertex(end,:) = vertex(1,:);
            obj.polygon = polyshape(vertex);
        end

        function vertex = flip(obj,vertex,vertElem,dimFig)
            nElem = obj.mesh.nelem;
            nodes = dimFig*nElem+1:vertElem*nElem;
            vertex(nodes,1) = - fliplr(vertex(1:dimFig*nElem,1)')';
            vertex(dimFig*nElem+1:vertElem*nElem,2) = fliplr(vertex(1:dimFig*nElem,2)')';
        end

        function plotFigure(obj)
            pgon = obj.polygon;
            figure(3)
            clf
            plot(pgon);
            axis([-1 1 0 1])
            grid on
            grid minor
            title('Clamped-Clamped Column Profile (2D)','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
        end

    end 
    
        methods (Access = protected)
        
        function initPlotting(obj)
            %obj.patchHandle = patch(obj.axes,'Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,...
              %'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
              set(obj.axes,'ALim',[0, 1],'XTick',[],'YTick',[]);
              % obj.plot();
              obj.BCplotter.plot();
        end

        function createCamera(obj)
            nullAxes = axes;
            obj.cam = Camera_Null(nullAxes);
        end


    end
    
    methods (Access = protected, Static)
        
        function color = getColor()
            color = [0 0 0];
        end
        
    end
    
    
    methods (Access = private)
            
    end
    
end