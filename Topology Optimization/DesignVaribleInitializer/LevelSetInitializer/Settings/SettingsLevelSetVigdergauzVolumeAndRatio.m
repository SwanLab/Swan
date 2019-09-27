classdef SettingsLevelSetVigdergauzVolumeAndRatio < SettingsLevelSetCreator
    
    properties (Access = public)
        volumeMicro
        superEllipseRatio
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetVigdergauzVolumeAndRatio(varargin)
            obj.loadParams('paramsLevelSetCreator_VigdergauzVolumeAndRatio.json')            
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end