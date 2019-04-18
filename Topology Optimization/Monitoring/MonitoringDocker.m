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
        
        function obj = MonitoringDocker(cParams)
            obj.shallDisplayParams    = cParams.showOptParams;
            obj.shallDisplayDesignVar = cParams.plotting;
            obj.createMonitors(cParams);
        end
        
        function refresh(obj,x,it,cost,constraint,convergenceVars,hasFinished,istep,nstep)
            obj.paramsMonitor.refresh(it,cost,constraint,convergenceVars,hasFinished,istep,nstep);
            obj.designVarMonitor.refresh(x);
        end
        
    end
    
    methods (Access = private)
        
        function createMonitors(obj,cParams)
            settings = cParams.settings;
            designVar = cParams.designVar;
            convergenceVars = cParams.convergenceVars;
            obj.paramsMonitor = ParamsMonitorFactory.create(obj.shallDisplayParams,settings,convergenceVars);
            obj.designVarMonitor = DesignVarMonitorFactory().create(obj.shallDisplayDesignVar,settings,designVar);
        end
        
    end
    
end