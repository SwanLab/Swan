classdef DesignVarMonitor_AreaColumn < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarName = 'Area Column'
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_AreaColumn(cParams)
            obj@DesignVarMonitor_Abstract(cParams);
        end
        
        function plot(obj)
            z = obj.designVar.getColumnArea; 
            coord = obj.mesh.coord;
            nelem = obj.mesh.nelem;
            nnod = nelem+1;
            scale = 0.3;
            dimFig   = 2;
            vertElem = 4;
            vertex = zeros(vertElem*nelem+1,dimFig);
            for iNod = 1:nnod-1
                vertex(2*iNod-1,:)   = [scale*z(iNod)/2 coord(iNod)];
                vertex(2*iNod,:) = [scale*z(iNod)/2 coord(iNod+1)];
            end
            vertex(dimFig*nelem+1:vertElem*nelem,1) = - fliplr(vertex(1:dimFig*nelem,1)')';
            vertex(dimFig*nelem+1:vertElem*nelem,2) = fliplr(vertex(1:dimFig*nelem,2)')';
            vertex(end,:) = vertex(1,:);
            pgon = polyshape(vertex);
            figure(1)
            clf
            plot(pgon);
            grid on
            grid minor
            title('Clamped-Clamped Column Profile (2D)','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            
        end
        
    end
    
        methods (Access = protected)
        
        function initPlotting(obj)
%            obj.patchHandle = patch(obj.axes,'Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,...
 %              'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
            set(obj.axes,'ALim',[0, 1],'XTick',[],'YTick',[]);
            % obj.plot();
           obj.BCplotter.plot();
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