classdef SettingsMeshUnfitted < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMeshUnfitted.json'
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
                case 2
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
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createBackgroundMesh();
            obj.createBackgroundInterpolation();
        end
        
        function createBackgroundMesh(obj)
            if ischar(obj.meshBackground)
                obj.meshBackground = Mesh_GiD(obj.meshBackground);
            end
        end
        
        function createBackgroundInterpolation(obj)
            obj.interpolationBackground = Interpolation.create(obj.meshBackground,'LINEAR');
        end
        
    end
    
end