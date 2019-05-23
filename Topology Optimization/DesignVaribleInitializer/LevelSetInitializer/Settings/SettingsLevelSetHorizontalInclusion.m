classdef SettingsLevelSetHorizontalInclusion < SettingsLevelSetCreator
    
    properties (Access = public)
        widthH
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetHorizontalInclusion(varargin)
            obj.loadParams('paramsLevelSetCreator_HorizontalInclusion.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end