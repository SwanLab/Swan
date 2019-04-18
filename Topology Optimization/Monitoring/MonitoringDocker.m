classdef MonitoringDocker < handle
    
    properties (Access = private)
        paramsMonitor
        designVarMonitor
    end
    
    methods (Access = public)
        
        function obj = MonitoringDocker(cParams)
            obj.createMonitors(cParams);
        end
        
        function refresh(obj,it,hasFinished,istep,nstep)
            obj.paramsMonitor.refresh(it,hasFinished,istep,nstep);
            obj.designVarMonitor.refresh();
        end
        
    end
    
    methods (Access = private)
        
        function createMonitors(obj,cParams)
            obj.paramsMonitor = ParamsMonitorFactory.create(cParams.settingsParamsMonitor);
            obj.designVarMonitor = DesignVarMonitorFactory().create(cParams.settingsDesignVarMonitor);
        end
        
    end
    
end