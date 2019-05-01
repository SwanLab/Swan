classdef SettingsOptimizer < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsOptimizer'
    end
    
    properties (Access = public)
        nconstr
        target_parameters
        constraint_case
        name
        maxiter
        
        designVar
        cost
        constraint
        
        shallPrint
        printMode
        
        uncOptimizerSettings
        settingsMonitor
    end
    
    methods (Access = public)
        
        function obj = SettingsOptimizer(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
        function setupSettingsMonitor(obj,settings,isOld)
            if isOld
                obj.settingsMonitor.showOptParams               = settings.monitoring;
                obj.settingsMonitor.refreshInterval             = settings.monitoring_interval;
                obj.settingsMonitor.shallDisplayDesignVar       = settings.plotting;
                obj.settingsMonitor.shallShowBoundaryConditions = settings.showBC;
            end
            
            obj.settingsMonitor.problemID                   = settings.case_file;
            obj.settingsMonitor.costFuncNames               = settings.cost;
            obj.settingsMonitor.costWeights                 = settings.weights;
            obj.settingsMonitor.constraintFuncs             = settings.constraint;
            obj.settingsMonitor.dim                         = settings.pdim;
        end
        
    end
    
end