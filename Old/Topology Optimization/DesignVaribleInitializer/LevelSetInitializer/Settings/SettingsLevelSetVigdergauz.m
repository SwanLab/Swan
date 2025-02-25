classdef SettingsLevelSetVigdergauz < SettingsLevelSetCreator
    
    properties (Access = public)
       vigdergauzDataBase
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetVigdergauz(varargin)
            obj.loadParams('paramsLevelSetCreator_Vigdergauz.json')            
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end