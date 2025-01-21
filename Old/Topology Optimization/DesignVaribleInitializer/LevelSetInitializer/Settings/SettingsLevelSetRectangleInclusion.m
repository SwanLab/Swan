classdef SettingsLevelSetRectangleInclusion < SettingsLevelSetCreator
    
    properties (Access = public)
        widthH
        widthV
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetRectangleInclusion(varargin)
            obj.loadParams('paramsLevelSetCreator_RectangleInclusion.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end