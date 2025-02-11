classdef Builder_DesignVarMonitor_2D < Builder_DesignVarMonitor_Abstract
    
    methods (Access = protected)
        
        function cam = createCamera(obj)
            cam = Camera_TopView(obj.monitor.axes);
        end
        
    end
    
end