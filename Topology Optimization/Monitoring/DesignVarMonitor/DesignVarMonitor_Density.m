classdef DesignVarMonitor_Density < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarName = 'Density - \rho';
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Density(mesh,showBC)
            obj@DesignVarMonitor_Abstract(mesh,showBC);
        end
        
        function plot(obj,rho)
            set(obj.patchHandle,'FaceVertexAlphaData',double(rho));
        end
        
    end
    
    methods (Access = protected)
        
        function initPlotting(obj)
            nnode = size(obj.mesh.coord,1);
            obj.patchHandle = patch(obj.axes,'Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,'FaceVertexAlphaData',zeros(nnode,1),...
                'FaceAlpha','flat','EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
            set(gca,'ALim',[0, 1],'XTick',[],'YTick',[]);
            
            obj.BCplotter.plot();
        end
        
    end
    
    methods (Access = protected, Static)
        
        function color = getColor()
            color = [0 0 0];
        end
        
    end
    
end