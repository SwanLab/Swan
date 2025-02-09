classdef SettingsMeshUnfitted < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMeshUnfitted.json'
    end
    
    properties (Access = public)
        backgroundMesh
        boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = SettingsMeshUnfitted(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
                case 2
                    obj.backgroundMesh = varargin{2};
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
        end
        
        function createBackgroundMesh(obj)
            if ischar(obj.backgroundMesh)
                fileName = obj.backgroundMesh;
                femReader = FemInputReaderGiD();
                s = femReader.read(fileName);
                obj.backgroundMesh = s.mesh;
            end
        end
        
    end
    
end