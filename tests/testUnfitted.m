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
            filename = fullfile('.','Input',obj.testName);
            settings = Settings(filename);
            settings.printChangingFilter = false;
            sett = SettingsTopOptProblem(settings);
            obj.topOpt = TopOpt_Problem(sett);
            obj.levelSet = obj.topOpt.designVariable.value;
        end
        
        function createMesh(obj)
            meshBackground = obj.topOpt.designVariable.mesh;
            interpolation = Interpolation.create(meshBackground,'LINEAR');
            settingsUnfitted = SettingsMeshUnfitted(obj.meshType,meshBackground,interpolation,obj.meshIncludeBoxContour);
            obj.mesh = Mesh_Unfitted_Factory.create(settingsUnfitted);
            obj.mesh.computeMesh(obj.levelSet);
        end
        
    end
    
end

