classdef DesignVarMonitor_Abstract < handle
    
    properties (Access = protected, Abstract)
        designVarName
    end
    
    properties (Access = protected)
        figHandle
        patchHandle
        designVar
        mesh
        bc
        cam
        BCplotter
    end
    
    properties (Access = private)
        dim
        showBC
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
        
        function obj = DesignVarMonitor_Abstract(cParams)
            obj.showBC    = cParams.showBC;
            obj.designVar = cParams.designVar;
            obj.mesh      = obj.designVar.mesh;
            obj.bc        = cParams.bc;
            obj.dim       = cParams.dim;
            obj.init();
        end
        
        function refresh(obj)
            obj.plot();
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
            obj.createCamera();
            obj.createBCplotter();
            obj.setupTheme();
            obj.initPlotting();
        end
        
        
        function initFrame(obj)
            obj.figHandle = figure;
            
            set(obj.figHandle,'Pointer','arrow','NumberTitle','off');
            title(obj.designVarName)
            
            hold on
            axis off
            axis equal
            
            obj.axes = obj.figHandle.Children;%CurrentAxes;
        end
        
        function setupTheme(obj)
            obj.figHandle.Color = 'white';
            colormap(obj.getColor());
        end
        
        function createCamera(obj)
            obj.cam = Camera_Null(obj.axes);
        end
        
        function createBCplotter(obj)
            obj.BCplotter = FactoryBoundayConditions.create(obj.showBC,obj.dim,obj.axes,obj.mesh,obj.bc);
        end
        
    end
    
end