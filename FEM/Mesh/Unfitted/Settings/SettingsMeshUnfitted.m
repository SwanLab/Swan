classdef SettingsMeshUnfitted < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMeshUnfitted'
    end
    
    properties (Access = public)
        unfittedType
        meshBackground
        interpolationBackground
    end
    
    methods (Access = public)
        
        function obj = SettingsMeshUnfitted(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
                case 3
                    obj.unfittedType = varargin{1};
                    obj.meshBackground = varargin{2};
                    obj.interpolationBackground = varargin{3};
            end
        end
        
    end
    
end