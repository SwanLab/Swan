classdef ParamsMonitorFactory <  handle
    
    methods (Access = public, Static)
        
        function monitor = create(shallDisplay,settings,convergenceVars)
            if shallDisplay
                monitor = ParamsMonitor(settings,convergenceVars);
            else
                monitor = ParamsMonitor_Null();
            end
        end
        
    end
    
end