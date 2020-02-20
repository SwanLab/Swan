classdef SettingsMeshUnfitted < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMeshUnfitted.json'
    end
    
    properties (Access = public)
        unfittedType
        meshBackground
        interpolationBackground
        includeBoxContour
        mesh
        type
        isInBoundary
    end
    
    methods (Access = public)
        
        function obj = SettingsMeshUnfitted(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
                case 2
                    obj.unfittedType = varargin{1};
                    mesh = varargin{2};
                    obj.meshBackground = mesh.innerMeshOLD;
                    obj.mesh           = mesh;
                case 3
                    disp('eis');
                case 4
                    disp('eis');
                    
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
                fileName = obj.meshBackground;
                femReader = FemInputReader_GiD();
                s = femReader.read(fileName);
                obj.meshBackground = s.mesh;
            end
        end
        
        function createBackgroundInterpolation(obj)
            inter = Interpolation.create(obj.meshBackground,'LINEAR');
            obj.interpolationBackground = inter;
        end
        
    end
    
end