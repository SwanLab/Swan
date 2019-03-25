classdef testUnfitted < test
    properties (Access = protected, Abstract)
        testName
        meshType
        meshIncludeBoxContour
    end
    
    properties (Access = protected)
        topOpt
        levelSet
        mesh
    end
    
    methods (Access = protected)
        function createTopOpt(obj)
            file_name_in = fullfile('.','Input',obj.testName);
            settings = Settings(file_name_in);
            settings.printChangingFilter = false;
            obj.topOpt = TopOpt_Problem(settings);
            obj.topOpt.preProcess;
            obj.levelSet = obj.topOpt.x;
        end
        
        function createMesh(obj)
            meshBackground = obj.topOpt.mesh;
            interpolation = Interpolation.create(meshBackground,'LINEAR');
            settingsUnfitted = SettingsMeshUnfitted(obj.meshType,meshBackground,interpolation);
            obj.mesh = Mesh_Unfitted_Factory.create(settingsUnfitted,'includeBoxContour',obj.meshIncludeBoxContour);
            obj.mesh.computeMesh(obj.levelSet);
        end
    end
end

