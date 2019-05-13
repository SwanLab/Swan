classdef SettingsIncrementalScheme < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsIncrementalScheme.json'
    end
    
    properties (Access = public)
        nSteps
        targetParamsSettings
        shallPrintIncremental
        mesh
        settings
    end
    
    methods (Access = public)
        
        function obj = SettingsIncrementalScheme(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createTargetParamsSettings();
        end
        
        function createTargetParamsSettings(obj)
            s = obj.targetParamsSettings;
            obj.targetParamsSettings = SettingsTargetParamsManager(s);
            if obj.settings.isOld
                obj.targetParamsSettings.VfracInitial = obj.settings.Vfrac_initial;
                obj.targetParamsSettings.VfracFinal = obj.settings.Vfrac_final;
                obj.targetParamsSettings.constrInitial = obj.settings.constr_initial;
                obj.targetParamsSettings.constrFinal = obj.settings.constr_final;
                obj.targetParamsSettings.optimalityInitial = obj.settings.optimality_initial;
                obj.targetParamsSettings.optimalityFinal = obj.settings.optimality_final;
            end
            obj.targetParamsSettings.epsilonInitial = obj.settings.epsilon_initial;
            obj.targetParamsSettings.epsilonFinal = obj.settings.epsilon_final;
            obj.targetParamsSettings.epsilonIsotropyInitial = obj.settings.epsilon_isotropy_initial;
            obj.targetParamsSettings.epsilonIsotropyFinal = obj.settings.epsilon_isotropy_final;
        end
        
    end
    
end