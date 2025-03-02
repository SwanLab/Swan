classdef SettingsLevelSetHorizontalFibers < SettingsLevelSetCreator
    
    properties
        levFib
        volume
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetHorizontalFibers(varargin)
            if nargin == 1
                obj.loadConfigFile(varargin{1})
            end
        end
        
    end
end