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
        
        femFileName
    end
    
    methods (Access = public)
        
        function obj = SettingsOptimizerUnconstrained(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
        function init(obj)
            obj.initScalarProductSettings();
            obj.initLineSearchSettings();
        end
        
    end
    
    methods (Access = private)
        
        function initScalarProductSettings(obj)
            obj.scalarProductSettings.filename = obj.femFileName;
        end
        
        function initLineSearchSettings(obj)
            s = obj.lineSearchSettings;
            obj.lineSearchSettings = SettingsLineSearch(s);
            s2.scalarProductSettings = obj.scalarProductSettings;
            s2.optimizerType  = obj.type;
            s2.filename       = obj.femFileName;
            obj.lineSearchSettings.loadParams(s2);
        end
        
    end
    
end