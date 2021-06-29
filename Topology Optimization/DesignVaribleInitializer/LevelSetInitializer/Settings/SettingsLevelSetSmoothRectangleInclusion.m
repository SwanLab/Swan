classdef SettingsLevelSetSmoothRectangleInclusion < SettingsLevelSetCreator
    
    properties (Access = public)
        widthH
        widthV
        pnorm
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetSmoothRectangleInclusion(varargin)
            obj.loadParams('paramsLevelSetCreator_SmoothRectangle.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end