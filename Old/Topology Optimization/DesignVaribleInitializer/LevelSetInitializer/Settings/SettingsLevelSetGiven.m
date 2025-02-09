classdef SettingsLevelSetGiven < SettingsLevelSetCreator
    
    properties (Access = public)
       % value
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetGiven(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end