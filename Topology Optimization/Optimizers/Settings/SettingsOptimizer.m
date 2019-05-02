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
        dualVariable
        cost
        constraint
        
        shallPrint
        printMode
        
        uncOptimizerSettings
        settingsMonitor
        incrementalScheme
        postProcessSettings
        historyPrinterSettings
    end
    
    methods (Access = public)
        
        function obj = SettingsOptimizer(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
        function setupSettingsHistoryPrinter(obj,fileName)
            obj.historyPrinterSettings.fileName   = fileName;
            obj.historyPrinterSettings.shallPrint = obj.shallPrint;
        end
        
        function setupSettingsPostProcess(obj)
            obj.postProcessSettings.shallPrint = obj.shallPrint;
            obj.postProcessSettings.printMode  = obj.printMode;
        end
        
    end
    
end