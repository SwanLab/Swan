classdef SettingsConstitutiveTensorFromVademecum < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsConstitutiveTensorFromVademecum'
    end
    
    properties (Access = public)
        fileName 
    end
    
    methods (Access = public)
        function obj = SettingsConstitutiveTensorFromVademecum(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
            end
        end
    end
    
end