classdef SettingsOptimizerUnconstrained < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsOptimizerUnconstrained.json'
    end
    
    properties (Access = public)
        type
        targetParameters
        designVariable
        lagrangian
        
        convergenceVars
        
        epsilon
        scalarProductSettings
        lineSearchSettings
        
        e2
        ub
        lb
        
        problemData
    end
    
    methods (Access = public)
        
        function obj = SettingsOptimizerUnconstrained(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
        function init(obj)
            obj.initProblemData();
            obj.initScalarProductSettings();
            obj.initLineSearchSettings();
        end
        
    end
    
    methods (Access = private)
        
        function initProblemData(obj)
            s = obj.problemData;
            obj.problemData = TopOptProblemDataContainer(s);
        end
        
        function initScalarProductSettings(obj)
            obj.scalarProductSettings.scale = obj.problemData.femData.scale;
            obj.scalarProductSettings.femSettings.fileName = obj.problemData.femData.fileName;
        end
        
        function initLineSearchSettings(obj)
            s = obj.lineSearchSettings;
            obj.lineSearchSettings = SettingsLineSearch(s);
            obj.lineSearchSettings.lineSearchInitiatorSettings.optimizerType  = obj.type;
            obj.lineSearchSettings.filename       = obj.problemData.femData.fileName;
        end
        
    end
    
end