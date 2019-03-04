classdef Builder_DesignVarMonitor_Abstract < handle
    
    properties (Access = protected)
        monitor
    end
    
    methods (Access = protected, Abstract)
        
        createCamera(obj)
        
    end
    
    methods (Access = public)
        
        function build(obj,monitor)
            obj.init(monitor)
            obj.setCamera();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,monitor)
            obj.monitor = monitor;
        end
        
        function setCamera(obj)
            cam = obj.createCamera();
            obj.monitor.setCamera(cam);
        end
        
    end
    
end