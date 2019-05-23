classdef SettingsLevelSetSeveralHoles < SettingsLevelSetCreator
    
    properties (Access = public)
        nHoles
        rHoles
        phaseHoles
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetSeveralHoles(varargin)
            obj.loadParams('paramsLevelSetCreator_SeveralHoles.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end