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
            mesh_background = obj.topOpt.mesh;
            obj.mesh = Mesh_Unfitted_Factory.create(obj.meshType,mesh_background,Interpolation.create(mesh_background,'LINEAR'),'includeBoxContour',obj.meshIncludeBoxContour);
            obj.mesh.computeMesh(obj.levelSet);
        end
    end
end

