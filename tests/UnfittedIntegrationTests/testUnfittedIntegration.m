classdef testUnfittedIntegration < testShowingError
    properties (Access = protected)
        tol = 6e-2;
        topOpt
        levelSet
        mesh
    end
    
    properties (Access = protected, Abstract)
        testName
        analiticalArea
        meshType
        meshIncludeBoxContour
    end
    
    properties (Access = private)
        areaAdim
    end
    
    methods (Access = protected)
        function obj = testUnfittedIntegration()
            obj.createTopOpt()
            obj.integrateSurface()
        end
        
        function createTopOpt(obj)
            file_name_in = strcat('./Input/',obj.testName);
            settings = Settings(file_name_in);
            settings.printChangingFilter = false;
            obj.topOpt = TopOpt_Problem(settings);
            obj.topOpt.preProcess;
        end
        
        function integrateSurface(obj)
            obj.levelSet = obj.topOpt.x;
            obj.createMesh();
            
            area = obj.computeGeometricalVariable();
            obj.areaAdim = area/obj.analiticalArea;
        end
        
        function createMesh(obj)
            mesh_background = obj.topOpt.mesh;
            obj.mesh = Mesh_Unfitted_Factory.create(obj.meshType,mesh_background,Interpolation.create(mesh_background,'LINEAR'),'includeBoxContour',obj.meshIncludeBoxContour);
            obj.mesh.computeMesh(obj.levelSet);
        end
        
        function M = computeGeometricalVariable(obj)
            M = obj.mesh.computeMass();
        end
        
        function computeError(obj)
            obj.error = abs(obj.areaAdim - 1);
        end
    end
end

