classdef DesignVarMonitor_Abstract < handle
    
    properties (Access = protected, Abstract)
        designVarName
    end
    
    properties (Access = protected)
        figHandle
        patchHandle
        mesh
        cam
    end
    
    properties (GetAccess = public, SetAccess = private)
        axes
    end
    
    methods (Access = public, Abstract)
        
        plot(obj)
        
    end
    
    methods (Access = protected, Abstract)
        
        initPlotting(obj)
        
    end
    
    methods (Access = protected, Static, Abstract)
        
        getColor()
        
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Abstract(mesh)
            obj.mesh = mesh;
            obj.init();
        end
        
        function refresh(obj,x)
            obj.plot(x);
            obj.cam.updateView();
            drawnow
        end
        
        function setCamera(obj,cam)
            obj.cam = cam;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.initFrame();
            obj.initPlotting();
            obj.setupTheme();
            obj.createCamera();
        end
        
        
        function initFrame(obj)
            obj.figHandle = figure;
            
            set(obj.figHandle,'Pointer','arrow','NumberTitle','off');
            title(obj.designVarName)
            
            axis off
            axis equal
            
            obj.axes = obj.figHandle.CurrentAxes;
        end
        
        function setupTheme(obj)
            obj.figHandle.Color = 'white';
            colormap(obj.getColor());
        end
        
        function createCamera(obj)
            obj.cam = Camera_Null(obj.axes);
        end
        
    end
    
end