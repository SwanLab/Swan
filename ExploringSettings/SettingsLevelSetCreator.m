classdef SettingsLevelSetCreator < AbstractSettings
    
    properties (Access = public)
        levelSetType
        ndim
        coord
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetCreator(varargin)
            obj.loadParams('paramsLevelSet');
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end