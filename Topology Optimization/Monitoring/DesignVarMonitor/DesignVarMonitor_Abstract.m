classdef DesignVarMonitor_Abstract < handle
    
    properties (Access = protected, Abstract)
        designVarName
    end
    
    properties (GetAccess = protected, SetAccess = private)
        figHandle
        patchHandle
        mesh
    end
    
    methods (Access = public, Abstract)
        
        plot(obj)
        
    end
    
    methods (Access = protected, Abstract)
        
        init(obj)
        
    end
    
    methods (Access = protected, Static, Abstract)
        
        getColor()
        
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Abstract(mesh)
            obj.mesh = mesh;
            obj.initFrame();
            obj.init();
        end
        
        function refresh(obj,x)
            obj.plot(x);
            drawnow
        end
        
    end
    
    methods (Access = private)
        
        function initFrame(obj)
            obj.figHandle = figure;
            set(obj.figHandle,'Pointer','arrow','NumberTitle','off');
            title(obj.designVarName)
            
            nnode = size(obj.mesh.coord,1);
            obj.patchHandle = patch('Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,'FaceVertexAlphaData',zeros(nnode,1),...
                'FaceAlpha','flat','EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
            set(gca,'ALim',[0, 1],'XTick',[],'YTick',[]);
            
            obj.setupTheme();
            axis off
            axis equal
        end
        
        function setupTheme(obj)
            obj.figHandle.Color = 'white';
            colormap(obj.getColor());
        end
        
    end
    
end