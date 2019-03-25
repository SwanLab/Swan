classdef SettingsMeshUnfitted < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMeshUnfitted'
    end
    
    properties (Access = public)
        unfittedType
        meshBackground
        interpolationBackground
        includeBoxContour
    end
    
    methods (Access = public)
        
        function obj = SettingsMeshUnfitted(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
                case 3
                    obj.unfittedType = varargin{1};
                    obj.meshBackground = varargin{2};
                    obj.interpolationBackground = Interpolation.create(obj.meshBackground,'LINEAR');
                case 3
                    obj.unfittedType = varargin{1};
                    obj.meshBackground = varargin{2};
                    obj.interpolationBackground = varargin{3};
                case 4
                    obj.unfittedType = varargin{1};
                    obj.meshBackground = varargin{2};
                    obj.interpolationBackground = varargin{3};
                    obj.includeBoxContour = varargin{4};
            end
        end
        
    end
    
end