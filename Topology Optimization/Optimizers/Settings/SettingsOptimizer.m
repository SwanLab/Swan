classdef SettingsOptimizer < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsOptimizer.json'
    end
    
    properties (Access = public)
        type
        problemData
        targetParameters
        constraintCase
        nConstr
        maxIter
        
        designVar
        dualVariable
        cost
        constraint
        
        shallPrint
        printMode
        
        uncOptimizerSettings
        monitoringDockerSettings
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
            if isfield(cParams,'monitoringDockerSettings')
                s = cParams.monitoringDockerSettings;
            else
                s = [];
            end
            obj.monitoringDockerSettings = SettingsMonitoringDocker(s);
            
            obj.monitoringDockerSettings.optimizerName   = obj.type;
            obj.monitoringDockerSettings.problemID       = obj.problemData.caseFileName;
            obj.monitoringDockerSettings.dim             = obj.problemData.pdim;
            
            obj.monitoringDockerSettings.costFuncNames   = obj.problemData.costFunctions;
            obj.monitoringDockerSettings.costWeights     = obj.problemData.costWeights;
            obj.monitoringDockerSettings.constraintFuncs = obj.problemData.constraintFunctions;
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