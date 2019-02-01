classdef MonitoringDocker < handle
    
    properties (GetAccess = public, SetAccess = private)
        shallDisplayParams
        shallDisplayDesignVar
    end
    
    properties (Access = private)
        paramsMonitor
        designVarMonitor
    end
    
    methods (Access = public)
        
        function obj = MonitoringDocker(shallDisplayParams,shallDisplayDesignVar,settings,mesh)
            obj.shallDisplayParams = shallDisplayParams;
            obj.shallDisplayDesignVar = shallDisplayDesignVar;
            obj.createMonitors(settings,mesh);
        end
        
        function refresh(obj,x,it,cost,constraint,convVars,hasFinished,istep,nstep)
            obj.paramsMonitor.refresh(it,cost,constraint,convVars,hasFinished,istep,nstep);
            obj.designVarMonitor.refresh(x);
        end
        
    end
    
    methods (Access = private)
        
        function createMonitors(obj,settings,mesh)
            obj.paramsMonitor = ParamsMonitorFactory.create(obj.shallDisplayParams,settings);
            obj.designVarMonitor = DesignVarMonitorFactory().create(obj.shallDisplayDesignVar,settings,mesh);
        end
        
    end
    
end