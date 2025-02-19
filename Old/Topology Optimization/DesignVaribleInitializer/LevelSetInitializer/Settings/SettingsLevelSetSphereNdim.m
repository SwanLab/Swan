classdef SettingsLevelSetSphereNdim < SettingsLevelSetCreator
    
    properties (Access = public)
        fracRadius
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetSphereNdim(varargin)
            obj.loadParams('paramsLevelSetCreator_SphereNdim.json')
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end