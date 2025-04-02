classdef SettingsLevelSetSmoothSquareInclusion < SettingsLevelSetCreator
    
    properties (Access = public)
        widthSquare
        pnorm
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetSmoothSquareInclusion(varargin)
            obj.loadParams('paramsLevelSetCreator_SmoothSquare.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end