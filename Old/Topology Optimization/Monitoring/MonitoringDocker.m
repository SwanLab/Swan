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
        
        function f = saveMonitorFigure(obj)
            f = obj.paramsMonitor.saveFigure();
        end
        
    end
    
    methods (Access = private)
        
        function createMonitors(obj,cParams)
            obj.createParamsMonitor(cParams);
            obj.createDesignVarMonitor(cParams);
        end
        
        function createParamsMonitor(obj,cParams)
            s.showOptParams         = cParams.showOptParams;
            s.refreshInterval       = cParams.refreshInterval;
            s.problemID             = cParams.problemID;
            s.costFuncNames         = cParams.costFuncNames;
            s.costWeights           = cParams.costWeights;
            s.constraintFuncs       = cParams.constraintFuncs;
            s.optimizerNames        = cParams.optimizerNames;
            s.designVariable      = cParams.designVariable;
            s.dualVariable        = cParams.dualVariable;
            s.cost                = cParams.cost;
            s.constraint          = cParams.constraint;
            s.convergenceVars     = cParams.convergenceVars;
            obj.paramsMonitor     = ParamsMonitorFactory.create(s);
        end
        
        function createDesignVarMonitor(obj,cParams)
            s.shallDisplay   = cParams.shallDisplayDesignVar;
            s.showBC         = cParams.shallShowBoundaryConditions;
            s.bc             = cParams.boundaryConditions;
            s.designVariable = cParams.designVariable;
            s.optimizerNames  = cParams.optimizerNames;
            s.dim            = cParams.dim;
            s.scale          = cParams.scale;
            s.mesh           = cParams.mesh;
            f = DesignVarMonitorFactory;
            obj.designVarMonitor = f.create(s);
        end
        
    end
    
end