classdef SettingsOptimizer < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsOptimizer.json'
    end
    
    properties (Access = public)
        problemData
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
        type
        
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
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
%             obj.initSettingsMonitorDocker();
        end
        
    end
    
    methods (Access = public)
        
        function initSettingsMonitorDocker(obj,cParams)
            s = cParams.settingsMonitor;
            
            s.optimizerName   = obj.name;
            s.problemID       = obj.problemData.caseFileName;
            s.dim             = obj.problemData.pdim;
            
            s.costFuncNames   = obj.problemData.costFunctions;
            s.costWeights     = obj.problemData.costWeights;
            s.constraintFuncs = obj.problemData.constraintFunctions;
            
            obj.settingsMonitor = SettingsMonitoringDocker(s);
        end
        
                
        function initSettingsHistoryPrinter(obj,fileName)
            obj.historyPrinterSettings.fileName   = fileName;
            obj.historyPrinterSettings.shallPrint = obj.shallPrint;
        end
        
        function initSettingsPostProcess(obj)
            obj.postProcessSettings.shallPrint = obj.shallPrint;
            obj.postProcessSettings.printMode  = obj.printMode;
        end
        
    end
    
end