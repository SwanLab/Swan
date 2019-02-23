classdef ParamsMonitorFactory <  handle
    
    methods (Access = public, Static)
        
        function monitor = create(shallDisplay,settings)
            if shallDisplay
                monitor = ParamsMonitor(settings);
            else
                monitor = ParamsMonitor_Null;
            end
        end
        
    end
    
end