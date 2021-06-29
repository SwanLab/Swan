classdef Builder_DesignVarMonitor_3D < Builder_DesignVarMonitor_Abstract
    
    methods (Access = protected)
        
        function cam = createCamera(obj)
            cam = Camera_Rotatory(obj.monitor.axes);
        end
        
    end
    
end