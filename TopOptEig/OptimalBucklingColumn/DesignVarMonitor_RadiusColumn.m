classdef DesignVarMonitor_RadiusColumn < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarName = 'Radius Column'
    end
    
    properties (Access = private)
        polygon
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_RadiusColumn(cParams)
            obj@DesignVarMonitor_Abstract(cParams);
        end
        
        function plot(obj)
            obj.nIter= obj.nIter+1;
            scale = 0.3;
            obj.createPolygon(scale);
            obj.plotFigure();
            obj.plotDesignVariable();
            %obj.create3Dplot();
        end
        
    end

    methods (Access = private)

        function createPolygon(obj,scl)
            z = obj.sectionVariables.computeArea;
            coord = obj.mesh.coord;
            nnod     = obj.mesh.nelem+1;
            dimFig   = 2;
            vertElem = 4;
            vertex = zeros(vertElem*obj.mesh.nelem+1,dimFig);
            for iNod = 1:nnod-1
                vertex(2*iNod-1,:)   = [z(iNod)/2 coord(iNod)];
                vertex(2*iNod,:) = [z(iNod)/2 coord(iNod+1)];
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
            yMax = max(obj.mesh.coord);
            figure(3)
            clf
            plot(pgon);
            axis([-10 10 0 yMax])
            grid on
            grid minor
            title('Column Profile (2D)','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
        end

        function plotDesignVariable(obj)
            y = obj.sectionVariables.getSingleValue();
            xMax = max(obj.mesh.coord);
            x = linspace(0,xMax,length(y));
            figure(4)
            clf
            plot(x,y);
            axis([0 xMax 0 2])
            grid on
            grid minor
            title('Design Variable','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('R(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
        end

        function create3Dplot(obj)
            nElem = obj.sectionVariables.mesh.nelem;
            nVar = obj.sectionVariables.designVariable.nDesignVar;
            s.designVariableValue = obj.sectionVariables.designVariable.value(1:nVar*nElem);
            s.coord = obj.sectionVariables.mesh.coord;
            s.type = 'cylinderBuckling'; 
            plt = Plot3DBucklingColumn(s);
            plt.compute();
            obj.frames{obj.nIter} = plt.frame;
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