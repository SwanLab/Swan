classdef SettingsLevelSetCircleInclusion < SettingsLevelSetCreator
    
    properties (Access = public)
        fracRadius
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetCircleInclusion(varargin)
            obj.loadParams('paramsLevelSetCreator_CircleInclusion.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end