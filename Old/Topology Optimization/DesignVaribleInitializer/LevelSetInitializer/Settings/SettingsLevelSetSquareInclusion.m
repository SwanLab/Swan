classdef SettingsLevelSetSquareInclusion < SettingsLevelSetCreator
    
    properties (Access = public)
        widthSquare
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetSquareInclusion(varargin)
            obj.loadParams('paramsLevelSetCreator_SquareInclusion.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end