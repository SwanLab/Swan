classdef MonitoringDocker < handle
    
    properties (Access = private)
        paramsMonitor
        designVarMonitor
    end
    
    methods (Access = public)
        
        function obj = MonitoringDocker(shallDisplayParams,shallDisplayDesignVar,settings,mesh)
            obj.createMonitors(shallDisplayParams,shallDisplayDesignVar,settings,mesh);
        end
        
        function refresh(obj,x,it,cost,constraint,convVars,hasFinished,istep,nstep)
            obj.paramsMonitor.refresh(it,cost,constraint,convVars,hasFinished,istep,nstep);
            obj.designVarMonitor.refresh(x);
        end
        
    end
    
    methods (Access = private)
        
        function createMonitors(obj,shallDisplayParams,shallDisplayDesignVar,settings,mesh)
            obj.paramsMonitor = ParamsMonitorFactory.create(shallDisplayParams,settings);
            obj.designVarMonitor = DesignVarMonitorFactory().create(shallDisplayDesignVar,settings,mesh);
        end
        
    end
    
end