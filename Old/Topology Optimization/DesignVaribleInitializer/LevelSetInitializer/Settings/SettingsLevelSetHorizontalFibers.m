classdef SettingsLevelSetHorizontalFibers < SettingsLevelSetCreator
    
    properties (Access = public)
        levelOfFibers        
        volume
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetHorizontalFibers(varargin)
            obj.loadParams('paramsLevelSetCreator_HorizontalFibers.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end