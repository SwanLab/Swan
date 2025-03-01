classdef ParamsMonitorFactory <  handle
    
    methods (Access = public, Static)
        
        function monitor = create(cParams)
            shallDisplay = cParams.showOptParams;
            if shallDisplay
                monitor = ParamsMonitor(cParams);
            else
                monitor = ParamsMonitor_Null();
            end
        end
        
    end
    
end